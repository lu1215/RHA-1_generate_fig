import yaml
import os

def update_yaml(path, updates):
    with open(path, 'r') as f:
        config = yaml.safe_load(f) or {}

    for key_path, value in updates.items():
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
        yaml.dump(config, f, default_flow_style=False)

# =====================
# 主程式區域
# =====================
scatter_groups = ["OG-02", "OG-06", "OG-10", "OG-14"]
run_config_path = "sRNAanalyst/example/config/run_config.yml"
plot_config_path = "sRNAanalyst/example/config/plot_config.yml"
stylesheet_path = "sRNAanalyst/example/config/stylesheet.yml"
cwd = os.getcwd()
print(cwd)

style_updates = {}
style_updates["General.style"] = "white"
update_yaml(stylesheet_path, style_updates)

for i in range(0, len(scatter_groups), 2):
    og1 = scatter_groups[i]
    og2 = scatter_groups[i + 1]
    print("scatter analysis: {} and {}".format(og1, og2))


    # run_config 設定
    run_updates = {
        "Metagene.run": False,
        "Scatter.run": True
    }

    if i == 0:
        run_updates["Scatter.data[0].name"] = "WAGO-1_IP"
        run_updates["Scatter.data[1].name"] = "WAGO-1_IP_rha-1_MUT"
    else:
        run_updates["Scatter.data[0].name"] = "CSR-1_IP"
        run_updates["Scatter.data[1].name"] = "CSR-1_IP_rha-1_MUT"

    update_yaml(run_config_path, run_updates)

    # plot_config 設定
    plot_updates = {
        "Metagene.plot": False,
        "Scatter.plot": True
    }

    if i == 0:
        plot_updates["Scatter.data[0].name"] = "WAGO-1_IP"
        plot_updates["Scatter.data[1].name"] = "WAGO-1_IP_rha-1_MUT"
    else:
        plot_updates["Scatter.data[0].name"] = "CSR-1_IP"
        plot_updates["Scatter.data[1].name"] = "CSR-1_IP_rha-1_MUT"

    # 第一次分析（CSR-1_target）
    plot_updates["Scatter.filter[0].name"] = "CSR-1_target"
    update_yaml(plot_config_path, plot_updates)


    os.chdir("sRNAanalyst/example")

    os.system("sh script/22G_analyze.sh")
    if i == 0:
        os.system("mv output/analyze/fig/Scatter_0.png ../../output/fig_f1.png")
        # os.system("mv output/analyze/fig/Scatter_0.svg ../../output/fig_f1.svg")
    else:
        os.system("mv output/analyze/fig/Scatter_0.png ../../output/fig_e1.png")
        # os.system("mv output/analyze/fig/Scatter_0.svg ../../output/fig_e1.svg")

    os.chdir(cwd)
    # 第二次分析（WAGO-1_target）
    plot_updates["Scatter.filter[0].name"] = "WAGO-1_target"
    update_yaml(plot_config_path, plot_updates)

    os.chdir("sRNAanalyst/example")
    os.system("sh script/22G_analyze.sh")
    if i == 0:
        os.system("mv output/analyze/fig/Scatter_0.png ../../output/fig_f2.png")
        # os.system("mv output/analyze/fig/Scatter_0.svg ../../output/fig_f2.svg")
    else:
        os.system("mv output/analyze/fig/Scatter_0.png ../../output/fig_e2.png")
        # os.system("mv output/analyze/fig/Scatter_0.svg ../../output/fig_e2.svg")

    os.chdir(cwd)
