import json
from pymatgen.core import Structure
from pymatgen.core import Element

with open("example_json/ABECAL.json", "r") as f:
    data = json.load(f)

print(data["geometry"].keys())
struct = data["geometry"]["relaxed_crystal"]
s = Structure.from_dict(struct)
s.remove_species([Element("H")])

data["geometry"]["wrong_crystal"] = s.as_dict()

with open("example_json/ABECAL_example.json", "w") as f:
    json.dump(data, f, indent=4)

# idx_list = []
# for i, site in s.sites:
#     if site.specie == Element("H"):
#         idx_list.append(i)


# import code
# code.interact(local=locals())