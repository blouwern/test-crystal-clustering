`dump_evt_info_with_time.py`: analyse to get any event with duplicate hits on a same modID and output to `*.log` files
`parse_time_diff.py`: parse text output of a .log file and draw statistical results as time_diffs.png
# Result shows that the duplicate hits have delta t at the nanosecond level
# This indicates that the duplicate hits brobably come from the same particle trajectory
# So, when read the raw file to get the edeps information, its better to reuse the `read_edep_of_each_evt.C` file to sum up the edeps, rather than filter the timeout hits in the `read_evt_edeps.py`
