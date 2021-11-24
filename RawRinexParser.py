import os
import re
from typing import List, Tuple, Union, Optional

OBSERVATION_FILE = os.path.join(os.getcwd(), "example_data", "940779338I_1.17O")
NAVIGATION_FILE = os.path.join(os.getcwd(), "example_data", "940779338I_1.17N")
END_OF_HEADER = 'END OF HEADER'
IONO_COEFFICIENT_LABELS = ['ION ALPHA', 'ION BETA']
TIME_OF_FIRST_OBS = 'TIME OF FIRST OBS'


def get_lines_from(file: str) -> Tuple[List[str], List[str]]:
    """
    Separate header from data in rinex files
    :param file: file to divide
    :return: separeted header and data
    """
    with open(file) as f:
        data_lines = f.read().split('\n')
        for i, line in enumerate(data_lines):
            if END_OF_HEADER in line:
                return data_lines[:i], data_lines[i + 1:]


def extract_numbers_from(line: str) -> List[Union[float]]:
    """
    Extract all numbers from one string using regex
    :param line: string to extract data
    :return: list of extracted numbers
    """
    match_number = re.compile('-? *[0-9]+\.?[0-9]*(?:[Ee] *[+-]? *[0-9]+)?')
    return [float(x) for x in re.findall(match_number, line)]


def extract_iono_coefficients(header_data: List[str]) -> List[float]:
    """
    Extract 8 coefficients allowing calculation of estimated ionospheric delays
    :param header_data: lines extracted from rinex navigation file
    :return: list of 8 coefficients
    """
    final_iono_list = []
    for line in header_data:
        if any([label in line for label in IONO_COEFFICIENT_LABELS]):
            [final_iono_list.append(item) for item in extract_numbers_from(line)]
    return final_iono_list


def extract_gps_time(header_data: List[str]) -> int:
    """
    Extract gps time from provided header data
    """
    for line in header_data:
        if TIME_OF_FIRST_OBS in line:
            nums = extract_numbers_from(line)
            return int(24 * 3600 + nums[3] * 3600 + nums[4] * 60 + nums[5])


def extract_distances_from(lines: List[str]) -> Optional[Tuple[List[Tuple[int, float]], str]]:
    """
    Extract pseudo distances to satellites based on provided
    lines. Return tuple with list of pseudo distances and string
    representing satellite header.
    """
    final_pseudo_distances = []
    satelites = [int(x) for x in extract_numbers_from(lines[0].split()[7])[1::]]
    start_id = 1
    for i in range(len(satelites)):
        nums = extract_numbers_from(lines[start_id])
        start_id += 2
        final_pseudo_distances.append((satelites[i], nums[0]))
    if len(final_pseudo_distances) < 1:
        return None
    else:
        return final_pseudo_distances, lines[0]


def find_pseudo_distance_blocks(observation_data: List[str]) -> List:
    """
    Return list of pseudo distance objects extracted from observation rinex file
    """
    last_checked_index = -1
    pseudo_distances = []
    for i, line in enumerate(observation_data):
        if i == len(observation_data) - 1:
            break
        if 'G' in line and last_checked_index < i:
            for j in range(i + 1, len(observation_data) - 1):
                if 'G' in observation_data[j]:
                    pseudos = extract_distances_from(observation_data[i:j - 1])
                    pseudo_distances.append(pseudos)
                    last_checked_index = j
                    break
    return pseudo_distances


def timer_parameters_of(number_of_satellites: int, navigation_data: List[str]) -> List[Tuple[int, List[float]]]:
    """
    Return list of tuples consist from satellite number and list of timer params
    """
    final_tuple_list = []
    actual_index = 0
    ix = 0
    while ix < len(navigation_data) - 1:
        indexes = [(8 * i + actual_index) for i in range(number_of_satellites)]
        actual_index += number_of_satellites * 8
        ix = actual_index
        for index in indexes:
            timer = extract_timer_param_starting_from(index, navigation_data)
            final_tuple_list.append(timer)
        return final_tuple_list


def extract_timer_param_starting_from(line: int, data: List[str]) -> Tuple[int, List[float]]:
    """
    Return list tuple consist from satellite number and list of timer params
    """
    satellite_num = int(extract_numbers_from(data[line])[0])
    timer_params = extract_numbers_from(data[line + 1])[1::]
    timer_params.extend(extract_numbers_from(data[line + 2]))
    timer_params.extend(extract_numbers_from(data[line + 3]))
    timer_params.extend(extract_numbers_from(data[line + 4]))
    timer_params.extend([extract_numbers_from(data[line + 5])[0], 0.0, 0.0])
    timer_params.extend(extract_numbers_from(data[line])[-3::])
    return satellite_num, timer_params


def main() -> None:
    nav_header, nav_data = get_lines_from(NAVIGATION_FILE)
    obs_header, obs_data = get_lines_from(OBSERVATION_FILE)
    pseudo_blocks = find_pseudo_distance_blocks(obs_data)
    satellites_count = len(pseudo_blocks[0][0])
    timer_blocks = timer_parameters_of(satellites_count, nav_data)
    time = extract_gps_time(obs_header)
    iono = extract_iono_coefficients(nav_header)
    pseudo_dist_header = pseudo_blocks[0][1]
    pseudo_distances = pseudo_blocks[0][0]
    timer_params = timer_blocks[0:4]

    with open("ouput_rinex.txt", "w") as f:
        f.write("GPS czas odbiornika\n")
        f.write(str(time) + "\n")
        f.write("Iono coefficients" + "\n")
        for io in iono:
            f.write(str(io) + "\n")
        f.write(f"Pseudoodleglosci {pseudo_dist_header}" + "\n")
        for pseudo_distance in pseudo_distances:
            f.write(f"{pseudo_distance[0]} {pseudo_distance[1]}" + "\n")
        f.write("0" + "\n")
        f.write("Efemerydy parametry zegarow kolejnych satelitow" + "\n")
        for timer_param in timer_params:
            f.write(str(timer_param[0]) + "\n")
            for tick in timer_param[1]:
                f.write(str(tick) + "\n")


if __name__ == '__main__':
    main()
