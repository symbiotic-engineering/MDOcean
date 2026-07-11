import json
import pathlib
import collections

key_list = collections.defaultdict(list)
end_files = pathlib.Path('.').rglob('end.json')
number_files = pathlib.Path('.').rglob('numbers.json')
files = list(number_files) + list(end_files)

print(f"Checking the following json files for duplicate keys...")
[print(f"{str(f)}") for f in files]

[key_list[k].append(f) for f in files for k in json.load(open(f))]
keys = sorted(key_list.items())

ignore = ['analysis_time','postpro_time'] # these keys are duplicated across every end.json file intentionally
repeated_keys = [k for k,v in keys if len(v)>1 and k not in ignore]

if repeated_keys:
    [ print(f'{k} in {len(v)} files:\n\t\t\t{ ", ".join([str(f) for f in v]) }') for k,v in keys if k in repeated_keys ]

    # keys that are ok to keep in numbers.json (not duplicated)
    numbers_non_repeated = [k for k,v in keys if len(v)==1 and 'numbers.json' in str(v[0])]
    print(f"\nKeys that appear only once, in numbers.json:\n{', '.join(numbers_non_repeated)}")

    raise ValueError("Duplicate json keys found")
else:
    print("Success: no duplicate json keys found.")