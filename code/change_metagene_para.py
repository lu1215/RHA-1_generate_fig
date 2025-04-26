import yaml
import os
import sys

def update_yaml(path, updates):
    # print("Loading YAML: " + path)
    with open(path, 'r') as f:
        config = yaml.safe_load(f) or {}

    # Apply updates
    for key_path, value in updates.items():
        # print("Updating {} = {}".format(key_path, value))
        keys = key_path.split('.')
        ref = config
        for k in keys[:-1]:
            if k.endswith(']'):
                k, idx = k[:-1].split('[')
                idx = int(idx)
                ref = ref[k][idx]
            else:
                ref = ref[k]

        last_key = keys[-1]
        if last_key.endswith(']'):
            last_key, idx = last_key[:-1].split('[')
            idx = int(idx)
            ref[last_key][idx] = value
        else:
            ref[last_key] = value

    with open(path, 'w') as f:
        # print("Writing back to: " + path)
        yaml.dump(config, f, default_flow_style=False)


# 模擬 shell 中的 $i 變數
i = int(sys.argv[1]) if len(sys.argv) > 1 else 0
current_dir = os.path.abspath("sRNAanalyst/example")

# 修改 run_config.yml
run_config_updates = {
    "Metagene.run": True,
    "Scatter.run": False,
    "BASE_DIR": current_dir,
}

if i == 0:
    run_config_updates.update({
        "Metagene.data[0].name": "WAGO-1_IP",
        "Metagene.data[1].name": "WAGO-1_IP_rha-1_MUT"
    })
else:
    run_config_updates.update({
        "Metagene.data[0].name": "CSR-1_IP",
        "Metagene.data[1].name": "CSR-1_IP_rha-1_MUT"
    })

update_yaml("sRNAanalyst/example/config/run_config.yml", run_config_updates)

# 修改 plot_config.yml
plot_config_updates = {
    "Metagene.plot": True,
    "Scatter.plot": False,
    "BASE_DIR": current_dir,
}

if i == 0:
    plot_config_updates.update({
        "Metagene.data[0].name": "WAGO-1_IP",
        "Metagene.data[1].name": "WAGO-1_IP_rha-1_MUT"
    })
else:
    plot_config_updates.update({
        "Metagene.data[0].name": "CSR-1_IP",
        "Metagene.data[1].name": "CSR-1_IP_rha-1_MUT"
    })

update_yaml("sRNAanalyst/example/config/plot_config.yml", plot_config_updates)

## change stylesheet.yml
# if i == 0:
# original 
# style: darkgrid
# color: coolwarm
stylesheet_updates = {
    "General.style": "whitegrid",
    # "General.color": "coolwarm"
}
update_yaml("sRNAanalyst/example/config/stylesheet.yml", stylesheet_updates)