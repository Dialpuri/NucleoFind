import logging


def calculate_n_from_max_memory(max_memory: float):
    """
    Calculates the number of processing cores (n) based on the provided maximum memory
    constraint and a heuristic mapping of memory per core.

    This function uses a predefined heuristic mapping of memory per core to determine
    the number of cores that can fit within the specified maximum memory. The input
    memory value is expected to be a string representing the maximum available memory.

    Parameters:
        max_memory (str): The maximum allowable memory in string format.

    Returns:
        int: The calculated number of cores based on the memory constraint.
    """
    if not max_memory:
        return 1

    try:
        max_memory = float(max_memory)
    except ValueError:
        logging.error(
            f"Invalid max_memory value: {max_memory}. Please provide a valid number in GB, ignoring for now."
        )
        return 1

    # Check to ensure the max_memory isn't too big (i.e. passed in MB rather than GB)
    if max_memory > 1024:
        max_memory = max_memory // 1024
        logging.warning(
            f"NucleoFind has converted the max_memory argument to GB. Max Memory = {max_memory} GB"
        )

    if max_memory < 1:
        logging.error(
            f"Invalid max_memory value: {max_memory}. Please provide a valid number in GB, ignoring for now."
        )
        return 1

    if max_memory < 4:
        return 1

    max_cpus = (max_memory - 1.21) // 2.45
    return max_cpus
