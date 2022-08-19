import json
import time

timestamp = time.strftime("%Y%m%d-%H%M%S")

configured_samples = []
for key in config.keys():
    if not key.startswith("sample_description"):
        continue
    sample = key.split("_", 2)[-1]
    configured_samples.append(sample)

if configured_samples:
    second_dump = "config_{}_{}.json".format(timestamp, "_".join(sorted(configured_samples)))
else:
    second_dump = "config_{}.json".format(timestamp)

with open(output[0], "w") as fake:
    _ = fake.write(second_dump + "\n(Full configuration dump)")

with open(second_dump, "w") as dump:
    json.dump(
        config,
        dump,
        ensure_ascii=True,
        indent=2,
        sort_keys=True,
    )
