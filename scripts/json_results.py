import numpy as np
import ROOT
import json
import logging

folder = "Vega_sWeights"
asym_dict = {}
asym_dict["block1"] = {}
# for w in range(1,8):
for n in [0,1,2,3,4,5,6,7,8,9,10,11,15,16]:
    label = f"dalitz_bin_{n}_hlt1_kp"
    asym_dict["block1"][label] = {}
    
    json_file = f"./{folder}/BIN{n}/a_raw.json"

    try:
        with open(json_file, 'r') as f:
            results = json.load(f)

    except FileNotFoundError:
        print(f"Error: {json_file} not found - skipping bin {n}")
        continue  # Skip this bin instead of exiting
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON format in {json_file}")
        print(f"JSON decode error: {e}")
        continue  # Skip this bin instead of exiting

    iteration = f"5w_bin{n}"
    
    if iteration in results:
        Araw = results[iteration].get("Araw", "")
        asym_dict["block1"][label] = Araw
        print(f"Found Araw for {iteration}: {Araw}")
    else:
        print(f"Warning: {iteration} not found in {json_file}")
        print(f"Available keys: {list(results.keys())}")
        asym_dict["block1"][label] = None
    
        
with open(f"./results.json", "w") as f:
    json_str = json.dumps(asym_dict)
    f.write(json_str)