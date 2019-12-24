import sys
import os
import re

from collections import OrderedDict, namedtuple, defaultdict

from pyteomics.auxiliary import unitfloat

import numpy as np


def _try_number(string):
    try:
        x = float(string)
        return x
    except (TypeError, ValueError):
        return string


_DEFAULT_DLL_PATH = os.path.join(
    os.path.dirname(
        os.path.realpath(__file__)),
    "_vendor",
    "ThermoRawFileReader_3_0_41",
    "Libraries")


# late binding imports

Business = None
_RawFileReader = None
clr = None
NullReferenceException = Exception
Marshal = None
IntPtr = None
Int64 = None


def is_thermo_raw_file(path):
    '''Detect whether or not the file referenced by ``path``
    is a Thermo RAW file.

    Parameters
    ----------
    path: :class:`str`
        The path to test

    Returns
    -------
    :class:`bool`:
        Whether or not the file is a Thermo RAW file.
    '''
    if not _test_dll_loaded():
        try:
            register_dll()
        except ImportError:
            return False
    try:
        source = _RawFileReader.RawFileReaderAdapter.FileFactory(path)
        source.SelectInstrument(Business.Device.MS, 1)
        return True
    except NullReferenceException:   # pylint: disable=broad-except
        return False


def determine_if_available():
    '''Checks whether or not the COM-based Thermo
    RAW file reading feature is available.

    This is done by attempting to instantiate the
    COM-provided object, which queries the Windows
    registry for the MSFileReader.dll.

    Returns
    -------
    :class:`bool`:
        Whether or not the feature is enabled.
    '''
    try:
        return _register_dll([_DEFAULT_DLL_PATH])
    except (OSError, ImportError):
        return False


def _register_dll(search_paths=None):
    '''Start the Common Language Runtime interop service by importing
    the :mod:`clr` module from Pythonnet, and then populate the global
    names referring to .NET entities, and finally attempt to locate the
    ThermoRawFileReader DLLs by searching alogn ``search_paths``.

    Parameters
    ----------
    search_paths: list
        The paths to check along for the ThermoRawFileReader DLL bundle.

    Returns
    -------
    :class:`bool`:
        Whether or not the .NET library successfully loaded
    '''
    if search_paths is None:
        search_paths = []
    global _RawFileReader, Business, clr, NullReferenceException   # pylint: disable=global-statement
    global Marshal, IntPtr, Int64   # pylint: disable=global-statement
    if _test_dll_loaded():
        return True
    try:
        import clr  # pylint: disable=redefined-outer-name
        from System import NullReferenceException  # pylint: disable=redefined-outer-name
        clr.AddReference("System.Runtime")
        clr.AddReference("System.Runtime.InteropServices")
        from System import IntPtr, Int64  # pylint: disable=redefined-outer-name
        from System.Runtime.InteropServices import Marshal  # pylint: disable=redefined-outer-name
    except ImportError:
        return False
    for path in search_paths:
        sys.path.append(path)
        try:
            clr.AddReference('ThermoFisher.CommonCore.RawFileReader')
            clr.AddReference('ThermoFisher.CommonCore.Data')
        except OSError:
            continue
        try:
            import ThermoFisher.CommonCore.Data.Business as Business  # pylint: disable=redefined-outer-name
            import ThermoFisher.CommonCore.RawFileReader as _RawFileReader  # pylint: disable=redefined-outer-name
        except ImportError:
            continue
    return _test_dll_loaded()


def register_dll(search_paths=None):
    '''Register the location of the Thermo RawFileReader DLL bundle with
    the Common Language Runtime interop system and load the .NET symbols
    used by this feature.

    Parameters
    ----------
    search_paths: list
        The paths to check along for the ThermoRawFileReader DLL bundle.

    '''
    if search_paths is None:
        search_paths = []
    search_paths = list(search_paths)
    search_paths.append(_DEFAULT_DLL_PATH)
    loaded = _register_dll(search_paths)
    if not loaded:
        msg = '''The ThermoFisher.CommonCore libraries could not be located and loaded.'''
        raise ImportError(msg)


def _test_dll_loaded():
    return _RawFileReader is not None


def _copy_double_array(src):
    '''A quick and dirty implementation of the fourth technique shown in
    https://mail.python.org/pipermail/pythondotnet/2014-May/001525.html for
    copying a .NET Array[Double] to a NumPy ndarray[np.float64] via a raw
    memory copy.

    ``int_ptr_tp`` must be an integer type that can hold a pointer. On Python 2
    this is :class:`long`, and on Python 3 it is :class:`int`.
    '''
    dest = np.empty(len(src), dtype=np.float64)
    Marshal.Copy(
        src, 0,
        IntPtr.__overloads__[Int64](dest.__array_interface__['data'][0]),
        len(src))
    return dest


IsolationWindow = namedtuple("IsolationWindow", ['lower', 'target', 'upper'])
PrecursorInformation = namedtuple("PrecursorInformation", [
                                  'mz', 'intensity', 'charge', 'precursor_scan_id',
                                  'source', 'product_scan_id'])
ActivationInformation = namedtuple("ActivationInformation", ['method', 'energy'])
ScanWindow = namedtuple("ScanWindow", ['lower', 'upper'])


class ScanEventInformation(namedtuple('ScanEventInformation', ['start_time', 'window_list', 'injection_time', 'traits'])):
    __slots__ = ()

    def pack(self):
        d = {
            "scan start time": unitfloat(self.start_time, 'minute'),
            "ion injection time": unitfloat(self.injection_time, 'millisecond'),
        }
        d.update(self.traits)
        window_list = [
            {
                'scan window lower limit': unitfloat(sw.lower, 'm/z'),
                'scan window upper limit': unitfloat(sw.upper, 'm/z')
            }
            for sw in self.window_list
        ]
        d['scanWindowList'] = {
            "count": len(window_list),
            "scanWindow": window_list
        }
        return d

ScanAcquisitionInformation = namedtuple(
    "ScanAcquisitionInformation", ['combination', 'scan_list'])


class RawReaderInterface(object):

    def _scan_arrays(self, scan):
        scan_number = scan.scan_number + 1
        stats = self._source.GetScanStatsForScanNumber(scan_number)
        segscan = self._source.GetSegmentedScanFromScanNumber(
            scan_number, stats)
        mzs = _copy_double_array(segscan.Positions)
        inten = _copy_double_array(segscan.Intensities)
        return mzs, inten

    def _pick_peaks_vendor(self, scan, *args, **kwargs):
        scan_info = Business.Scan.FromFile(self._source, scan.scan_number + 1)
        if scan_info.HasCentroidStream:
            stream = self._source.GetCentroidStream(scan.scan_number + 1, 0)
            mzs = stream.Masses
            intens = stream.Intensities
            mzs = _copy_double_array(mzs)
            intens = _copy_double_array(intens)

            return mzs, intens
        else:
            raise NotImplementedError()

    def _scan_id(self, scan):
        scan_number = scan.scan_number
        return _make_id(scan_number + 1)

    def _is_profile(self, scan):
        return not self._source.IsCentroidScanFromScanNumber(
            scan.scan_number + 1)

    def _polarity(self, scan):
        filter_string = self._filter_string(scan)
        return filter_string.data['polarity']

    def _scan_title(self, scan):
        return "%s %r" % (self._scan_id(scan), self._filter_string(scan))

    def _filter_string(self, scan):
        scan_number = scan.scan_number
        return FilterString(self._source.GetFilterForScanNumber(scan_number + 1).Filter)

    def _scan_index(self, scan):
        scan_number = scan.scan_number
        return scan_number

    def _scan_time(self, scan):
        scan_number = scan.scan_number
        return self._source.RetentionTimeFromScanNumber(scan_number + 1)

    def _ms_level(self, scan):
        scan_number = scan.scan_number
        f = self._source.GetFilterForScanNumber(scan_number + 1)
        return f.MSOrder

    def _isolation_window(self, scan):
        scan_number = scan.scan_number
        filt = self._source.GetFilterForScanNumber(scan_number + 1)
        seq_index = filt.MSOrder - 2
        width = filt.GetIsolationWidth(seq_index)
        offset = filt.GetIsolationWidthOffset(seq_index)
        precursor_mz = filt.GetMass(seq_index)
        return IsolationWindow(width, precursor_mz + offset, width)

    def _trailer_values(self, scan):
        scan_number = scan.scan_number
        trailers = self._source.GetTrailerExtraInformation(scan_number + 1)
        return OrderedDict(zip([label.strip(":") for label in trailers.Labels], map(_try_number, trailers.Values)))

    def _infer_precursor_scan_number(self, scan):
        precursor_scan_number = None
        last_index = self._scan_index(scan) - 1
        current_level = self._ms_level(scan)
        i = 0
        while last_index >= 0 and i < 100:
            prev_scan = self.get_by_index(last_index)
            if prev_scan.ms_level >= current_level:
                last_index -= 1
            else:
                precursor_scan_number = prev_scan._data.scan_number
                break
            i += 1
        return precursor_scan_number

    def _precursor_information(self, scan):
        scan_number = scan.scan_number
        filt = self._source.GetFilterForScanNumber(scan_number + 1)
        precursor_mz = filt.GetMass(filt.MSOrder - 2)
        trailers = self._trailer_values(scan)
        _precursor_mz = float(trailers.get("Monoisotopic M/Z", 0))
        if _precursor_mz > 0:
            precursor_mz = _precursor_mz

        # imitate proteowizard's firmware bug correction
        isolation_window = self._isolation_window(scan)
        if (isolation_window.upper + isolation_window.lower) / 2 <= 2.0:
            if (isolation_window.target - 3.0 > precursor_mz) or (isolation_window.target + 2.5 < precursor_mz):
                precursor_mz = isolation_window.target
        elif precursor_mz not in isolation_window:
            precursor_mz = isolation_window.target
        charge = int(trailers.get("Charge State", 0))
        if charge == 0:
            charge = None
        inten = 0
        precursor_scan_number = None
        if precursor_scan_number is None:
            last_index = self._scan_index(scan) - 1
            current_level = self._ms_level(scan)
            i = 0
            while last_index >= 0 and i < 100:
                prev_scan = self.get_by_index(last_index)
                if prev_scan.ms_level >= current_level:
                    last_index -= 1
                else:
                    precursor_scan_number = prev_scan._data.scan_number
                    break
                i += 1
        if precursor_scan_number is not None:
            precursor_scan_id = self.get_by_index(
                precursor_scan_number).id
        return PrecursorInformation(
            precursor_mz, inten, charge, precursor_scan_id,
            source=self, product_scan_id=self._scan_id(scan))

    def _get_scan_segment(self, scan):
        trailer = self._trailer_values(scan)
        try:
            return int(trailer['Scan Segment'])
        except KeyError:
            return 1

    def _get_scan_event(self, scan):
        trailer = self._trailer_values(scan)
        try:
            return int(trailer['Scan Event'])
        except KeyError:
            return 1

    def _activation(self, scan):
        filter_string = self._filter_string(scan)
        tandem_sequence = filter_string.get("tandem_sequence")
        # If the tandem sequence exists, the last entry is the most recent tandem acquisition.
        # It will list contain one or more activation types. Alternatively, multiple activations
        # of the same precursor may exist in the list as separate events in the tandem sequence.
        if tandem_sequence is not None:
            activation_event = tandem_sequence[-1]
            activation_type = list(activation_event.get("activation_type"))
            has_supplemental_activation = filter_string.get(
                "supplemental_activation")

            if activation_type is not None:
                energy = list(activation_event.get("activation_energy"))
                if len(tandem_sequence) > 1:
                    prev_event = tandem_sequence[-2]
                    # Merge previous tandem sequences of the same precursor
                    if abs(prev_event['isolation_mz'] - activation_event['isolation_mz']) < 1e-3:
                        activation_type = list(prev_event.get(
                            "activation_type")) + activation_type
                        energy = list(prev_event.get(
                            "activation_energy")) + energy
                        has_supplemental_activation = True

                # if has_supplemental_activation and len(activation_type) > 1:
                #     activation_type.append(supplemental_term_map[
                #         dissociation_methods_map[activation_type[-1]]])
                if len(activation_type) == 1:
                    return ActivationInformation(activation_type[0], energy[0])
                else:
                    return ActivationInformation(activation_type, energy)
        return None

    def _acquisition_information(self, scan):
        fline = self._filter_string(scan)
        event = self._get_scan_event(scan)
        trailer_extras = self._trailer_values(scan)
        traits = {
            'preset scan configuration': event,
            'filter string': fline,
        }
        event = ScanEventInformation(
            self._scan_time(scan),
            injection_time=unitfloat(trailer_extras.get(
                'Ion Injection Time (ms)', 0.0), 'millisecond'),
            window_list=[ScanWindow(
                fline.get("scan_window")[0], fline.get("scan_window")[1])], traits=traits)
        return ScanAcquisitionInformation("no combination", [event])

    def _instrument_configuration(self, scan):
        fline = self._filter_string(scan)
        try:
            confid = self._analyzer_to_configuration_index[analyzer_map[fline.data.get(
                "analyzer")]]
            return self._instrument_config[confid]
        except KeyError:
            return None

    def _annotations(self, scan):
        fline = self._filter_string(scan)
        trailer_extras = self._trailer_values(scan)
        annots = {
            "filter string": fline,
        }
        microscans = trailer_extras.get("Micro Scan Count")
        if microscans is not None:
            annots['[Thermo Trailer Extra]Micro Scan Count'] = float(
                microscans)
        scan_segment = trailer_extras.get("Scan Segment")
        if scan_segment is not None:
            annots['[Thermo Trailer Extra]Scan Segment'] = int(scan_segment)
        scan_event = trailer_extras.get("Scan Event")
        if scan_event is not None:
            annots['[Thermo Trailer Extra]Scan Event'] = int(scan_event)
        mono_mz = float(trailer_extras.get("Monoisotopic M/Z", 0))
        if mono_mz is not None and mono_mz > 0:
            annots['[Thermo Trailer Extra]Monoisotopic M/Z'] = mono_mz
        hcd_ev = trailer_extras.get('HCD Energy eV')
        if hcd_ev is not None and hcd_ev > 0:
            annots['[Thermo Trailer Extra]HCD Energy eV'] = hcd_ev
        hcd_energies = trailer_extras.get('HCD Energy')
        if hcd_energies is not None and hcd_energies:
            annots['[Thermo Trailer Extra]HCD Energy'] = hcd_energies
        return annots


analyzer_pat = re.compile(
    r"(?P<mass_analyzer_type>ITMS|TQMS|SQMS|TOFMS|FTMS|SECTOR)")
polarity_pat = re.compile(r"(?P<polarity>[\+\-])")
point_type_pat = re.compile(r"(?P<point_type>[CP])")
ionization_pat = re.compile(
    r"(?P<ionization_type>EI|CI|FAB|APCI|ESI|APCI|NSI|TSP|FD|MALDI|GD)")
scan_type_pat = re.compile(r"(?P<scan_type>FULL|SIM|SRM|CRM|Z|Q1MS|Q3MS)")
ms_level_pat = re.compile(r" ms(?P<level>\d*) ")
activation_pat = re.compile(
    r"""(?:(?P<isolation_mz>\d+\.\d*)@
        (?P<activation_type>[a-z]+)
        (?P<activation_energy>\d*\.?\d*))""", re.VERBOSE)
activation_mode_pat = re.compile(
    r"""(?P<activation_type>[a-z]+)
        (?P<activation_energy>\d*\.\d*)""", re.VERBOSE)
scan_window_pat = re.compile(
    r"""
    \[(?P<scan_start>[0-9\.]+)-(?P<scan_end>[0-9\.]+)\]
    """, re.VERBOSE)

analyzer_map = {
    'FTMS': ("orbitrap"),
    "ITMS": ("ion trap"),
    "SQMS": ("quadrupole"),
    "TQMS": ("quadrupole"),
    "TOFMS": ("time-of-flight"),
    "SECTOR": ("magnetic sector")
}


ionization_map = {
    "EI": ("electron ionization"),
    "CI": ("chemical ionization"),
    "FAB": ("fast atom bombardment ionization"),
    "ESI": ("electrospray ionization"),
    "NSI": ("nanoelectrospray"),
    "APCI": ("atmospheric pressure chemical ionization"),
    "TSP": ("thermospray ionization"),
    "FD": ("field desorption"),
    "MALDI": ("matrix assisted laser desorption ionization"),
    "GD": ("glow discharge ionization"),
}


inlet_map = {
    "FAB": ("continuous flow fast atom bombardment"),
    "ESI": ("electrospray inlet"),
    "NSI": ("nanospray inlet"),
    "TSP": ("thermospray inlet"),
}


class ThermoRawScanPtr(object):
    def __init__(self, scan_number):
        self.scan_number = scan_number
        self.filter_string = None

    def validate(self, source):
        try:
            source._scan_time(self)
            return True
        except IOError:
            return False

    def __repr__(self):
        return "{self.__class__.__name__}({self.scan_number})".format(self=self)

    def pack(self, source):
        wrapper = ScanWrapper()
        wrapper._data = self
        wrapper['m/z array'], wrapper['intensity array'] = source._scan_arrays(self)
        wrapper['id'] = source._scan_id(self)
        wrapper['index'] = self.scan_number
        wrapper['ms level'] = source._ms_level(self)
        scan_acquisition = source._acquisition_information(self)
        wrapper['scanList'] = {
            'count': len(scan_acquisition.scan_list),
            scan_acquisition.combination: '',
            'scan': [se.pack() for se in scan_acquisition.scan_list]
        }
        return wrapper


class FilterString(str):
    def __init__(self, value):
        self.data = self._parse()

    def get(self, key):
        return self.data.get(key)

    def _parse(self):
        return filter_string_parser(self)


def filter_string_parser(line):
    """Parses instrument information from Thermo's filter string

    Parameters
    ----------
    line : str
        The filter string associated with a scan

    Returns
    -------
    dict
        Fields extracted from the filter string
    """
    words = line.upper().split(" ")
    values = dict()
    i = 0
    values['supplemental_activation'] = " sa " in line
    ms_level_info = ms_level_pat.search(line)
    if ms_level_info is not None:
        ms_level_data = ms_level_info.groupdict()
        level = ms_level_data.get("level")
        if level != "":
            parts = line[ms_level_info.end():].split(" ")
            tandem_sequence = []
            for part in parts:
                activation_info = activation_pat.search(part)
                if activation_info is not None:
                    activation_info = activation_info.groupdict()
                    activation_event = dict()
                    activation_event["isolation_mz"] = float(
                        activation_info['isolation_mz'])
                    activation_event["activation_type"] = [
                        activation_info['activation_type']]
                    activation_event["activation_energy"] = [
                        float(activation_info['activation_energy'])]
                    if part.count("@") > 1:
                        act_events = activation_mode_pat.finditer(part)
                        # discard the first match which we already recorded
                        next(act_events)
                        for match in act_events:
                            act_type, act_energy = match.groups()
                            act_energy = float(act_energy)
                            activation_event["activation_type"].append(
                                act_type)
                            activation_event['activation_energy'].append(
                                act_energy)
                    tandem_sequence.append(activation_event)
            values['ms_level'] = int(level)
            values['tandem_sequence'] = tandem_sequence

    scan_window_info = scan_window_pat.search(line)
    if scan_window_info is not None:
        values['scan_window'] = (
            float(scan_window_info.group(1)), float(scan_window_info.group(2)))

    try:
        word = words[i]
        i += 1
        analyzer_info = analyzer_pat.search(word)
        if analyzer_info is not None:
            values['analyzer'] = analyzer_info.group(0)
            word = words[i]
            i += 1
        polarity_info = polarity_pat.search(word)
        if polarity_info is not None:
            polarity_sigil = polarity_info.group(0)
            if polarity_sigil == "+":
                polarity = 1
            elif polarity_sigil == "-":
                polarity = -1
            else:
                polarity = 0
            values["polarity"] = polarity
            word = words[i]
            i += 1
        if word in "PC":
            if word == 'P':
                values['peak_mode'] = 'profile'
            else:
                values['peak_mode'] = 'centroid'
            word = words[i]
            i += 1
        ionization_info = ionization_pat.search(word)
        if ionization_info is not None:
            values['ionization'] = ionization_info.group(0)
            word = words[i]
            i += 1

        return values
    except IndexError:
        return values


_id_template = "controllerType=0 controllerNumber=1 scan="


def _make_id(scan_number):
    try:
        return "%s%d" % (_id_template, (scan_number))
    except TypeError:
        return None


def _parse_id(scan_id):
    return int(scan_id.replace(_id_template, ""))


class ScanWrapper(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    @property
    def scan_time(self):
        return 0

    @property
    def ms_level(self):
        return self['ms level']

    @property
    def index(self):
        return self['index']

    @property
    def scan_number(self):
        return self._data.scan_number


class ThermoRaw(RawReaderInterface):
    def __init__(self, source_file):
        if not _test_dll_loaded():
            register_dll()
        self._source = _RawFileReader.RawFileReaderAdapter.FileFactory(
            source_file)
        self._source.SelectInstrument(Business.Device.MS, 1)
        self.source_file = source_file
        # self._first_scan_time = self.get_by_index(0).scan_time
        # self._last_scan_time = self.get_by_index(self._source.RunHeaderEx.LastSpectrum - 1).scan_time
        self._producer = None
        self._scan_type_index = dict()
        self.make_iterator()
        self._index = self._pack_index()
        self._method = self._parse_method()
        self._get_instrument_info()

    def _build_scan_type_index(self):
        self.make_iterator()
        index = defaultdict(int)
        analyzer_counter = 1
        analyzer_confs = dict()
        for scan in self:  # pylint: disable=not-an-iterable
            index[scan.ms_level] += 1
            fline = self._filter_string(scan._data)
            analyzer = analyzer_map[fline.data['analyzer']]
            try:
                analyzer_confs[analyzer]
            except KeyError:
                analyzer_confs[analyzer] = analyzer_counter
                analyzer_counter += 1
        self.reset()
        self._scan_type_index = index
        self._analyzer_to_configuration_index = analyzer_confs

    def _get_instrument_info(self):
        scan = self.get_by_index(0)
        filter_string = self._filter_string(scan._data)
        ionization_label = filter_string.data.get("ionization")
        try:
            ionization = ionization_map[ionization_label]
        except KeyError:
            ionization = ionization_map['ESI']
        try:
            inlet = inlet_map[ionization_label]
        except KeyError:
            inlet = None

        source_group = set()
        source_group.add(ionization)
        if inlet is not None:
            source_group.add(inlet)
        configs = [source_group]
        return configs

    def make_iterator(self):
        self._producer = self._make_default_iterator()
        return self

    def __iter__(self):
        return iter(self._producer)

    def _make_scan_index_producer(self, start_index=None, start_time=None):
        if start_index is not None:
            return range(start_index, self._source.RunHeaderEx.LastSpectrum - 1)
        elif start_time is not None:
            start_index = self._scan_time_to_scan_number(start_time)
            while start_index != 0:
                scan = self.get_by_index(start_index)
                if scan.ms_level > 1:
                    start_index -= 1
                else:
                    break
            return range(start_index, self._source.RunHeaderEx.LastSpectrum - 1)
        else:
            return range(0, self._source.RunHeaderEx.LastSpectrum - 1)

    def _make_pointer_iterator(self, start_index=None, start_time=None):
        iterator = self._make_scan_index_producer(start_index, start_time)
        for i in iterator:
            yield ThermoRawScanPtr(i).pack(self)

    def _make_default_iterator(self):
        return self._make_pointer_iterator()

    @property
    def index(self):
        '''Accesses the scan index

        Returns
        -------
        :class:`collections.OrderedDict`
            Maps scan ID to index
        '''
        return self._index

    def __len__(self):
        return len(self.index)

    def __repr__(self):
        return "ThermoRaw(%r)" % (self.source_file)

    def close(self):
        '''Close the underlying file reader.
        '''
        if self._source is not None:
            self._source.Close()
            self._source = None

    def __del__(self):
        self.close()

    def reset(self):
        self.make_iterator()

    def _pack_index(self):
        index = OrderedDict()
        for sn in range(self._source.RunHeaderEx.FirstSpectrum - 1,
                        self._source.RunHeaderEx.LastSpectrum):
            index[_make_id(sn)] = sn
        return index

    def _get_instrument_model_name(self):
        return self._source.GetInstrumentData().Model

    def _get_instrument_serial_number(self):
        return self._source.GetInstrumentData().SerialNumber

    def _parse_method(self):
        try:
            method_count = self._source.InstrumentMethodsCount
            if method_count == 0:
                return ('')
            # the data acquisition method should be the last method
            return (self._source.GetInstrumentMethod(method_count - 1))
        except NullReferenceException:   # pylint: disable=broad-except
            return ('')

    def _scan_time_to_scan_number(self, scan_time):
        scan_number = self._source.ScanNumberFromRetentionTime(scan_time) - 1
        return scan_number

    def get_by_index(self, index):
        """Retrieve the scan object for the specified scan index.

        This internally calls :meth:`get_by_id` which will
        use its cache.

        Parameters
        ----------
        index: int
            The index to get the scan for

        Returns
        -------
        Scan
        """
        scan_number = int(index)
        package = ThermoRawScanPtr(scan_number)
        if not package.validate(self):
            raise IndexError(index)
        scan = package.pack(self)
        return scan

    def get_by_id(self, scan_id):
        """Retrieve the scan object for the specified scan id.

        If the scan object is still bound and in memory somewhere,
        a reference to that same object will be returned. Otherwise,
        a new object will be created.

        Parameters
        ----------
        scan_id : str
            The unique scan id value to be retrieved

        Returns
        -------
        Scan
        """
        scan_number = int(str(scan_id).replace(_id_template, '')) - 1
        package = ThermoRawScanPtr(scan_number)
        if not package.validate(self):
            raise KeyError(str(scan_id))
        scan = package.pack(self)
        return scan

    def get_by_time(self, time):
        """Retrieve the scan object for the specified scan time.

        This internally calls :meth:`get_by_id` which will
        use its cache.

        Parameters
        ----------
        time : float
            The time to get the nearest scan from

        Returns
        -------
        Scan
        """
        if time < self._first_scan_time:
            time = self._first_scan_time
        elif time > self._last_scan_time:
            time = self._last_scan_time
        scan_number = self._scan_time_to_scan_number(time)
        package = ThermoRawScanPtr(scan_number)
        scan = package.pack(self)
        return scan
