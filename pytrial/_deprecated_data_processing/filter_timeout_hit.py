from collections import defaultdict
import statistics


def filter_timeout_hit(evt_info_dict):
    """
    Args:
            evt_info_dict: {'ModID': [...], 'Edep': [...], 't': [...]}
            All 3 lists are aligned
    Returns:
            filterd dict: no ModID-duplicate timeout hit
    Logic:
     1. Identify groups with duplicate ModIDs, and exclude them all when calculating the benchmark time.
     2. Use the average time of the remaining unique ModIDs as the time reference t_ref.
     3. For each repeating group, retain the entry with the closest time to t_ref, and discard the rest.
    """
    # check items
    mod_ids = evt_info_dict.get("ModID", [])
    edeps = evt_info_dict.get("Edep", [])
    times = evt_info_dict.get("t", [])

    if not mod_ids:
        return {"ModID": [], "Edep": [], "t": []}

    # group<dict> structure: {modID:[<index_list>]}
    groups = defaultdict(list)
    for idx, mid in enumerate(mod_ids):
        groups[mid].append(idx)

    unique_indices = []
    duplicate_groups = []

    for mid, indices in groups.items():
        if len(indices) == 1:
            unique_indices.append(indices[0])
        else:
            duplicate_groups.append(indices)

    if unique_indices:

        t_ref = statistics.mean(times[i] for i in unique_indices)
    else:

        t_ref = statistics.mean(times)

    final_indices = list(unique_indices)

    for group in duplicate_groups:

        best_idx = min(group, key=lambda i: abs(times[i] - t_ref))
        final_indices.append(best_idx)

    final_indices.sort()

    return {
        "ModID": [mod_ids[i] for i in final_indices],
        "Edep": [edeps[i] for i in final_indices],
        "t": [times[i] for i in final_indices],
    }
