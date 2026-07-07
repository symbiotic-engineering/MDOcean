import json
import yaml
import glob

json_files = glob.glob('results/*/end.json')

all_data = {}

# Read from JSON file and write to YAML file
for json_file in json_files:
     with open(json_file, 'r') as f:
        data = json.load(f)
        for key, value in data.items():
            if key not in all_data:
                all_data[key] = value

with open('pubs/defense/quarto/_variables.yml', 'w') as yaml_file:
    yaml.dump(all_data, yaml_file, default_flow_style=False)