# # import sys
# # import pandas as pd
# # import matplotlib.pyplot as plt

# # def main(csv_file, transcript_name, transcript_length):
# #     # 讀取 CSV
# #     df = pd.read_csv(csv_file, comment='#')

# #     # 過濾指定 transcript
# #     df = df[df['ref_id'] == transcript_name]

# #     if df.empty:
# #         print("找不到指定的 transcript")
# #         return

# #     # 依照 init_pos 排序
# #     df = df.sort_values(by='init_pos')

# #     # 計算長條圖所需資料
# #     x = df['init_pos']
# #     widths = df['end_pos'] - df['init_pos']
# #     heights = df['read_count']

# #     # 畫長條圖
# #     plt.figure(figsize=(12, 6))
# #     plt.bar(x, heights, width=widths, align='edge', edgecolor='black')

# #     # 找最大 read_count 並標註
# #     max_row = df.loc[df['read_count'].idxmax()]
# #     max_x = (max_row['init_pos'] + max_row['end_pos']) // 2
# #     max_y = max_row['read_count']
# #     plt.text(max_x, max_y, str(max_y), fontsize=10, ha='center', va='bottom', color='red')

# #     # 軸與標題設定
# #     plt.xlabel('Transcript Position')
# #     plt.ylabel('Read Count')
# #     plt.title('{} mapped to PAR-6 gene({})'.format(csv_file.replace(".csv", '').split('/')[-1], transcript_name))
# #     plt.xlim(0, int(transcript_length))
# #     plt.tight_layout()

# #     # 儲存圖檔
# #     output_filename = "output/{}_22G_non_overlap.png".format(csv_file.replace(".csv", '').split('/')[-1])
# #     plt.savefig(output_filename)
# #     print("save graph to:", output_filename)

# #     # 顯示圖
# #     plt.show()

# # if __name__ == '__main__':
# #     if len(sys.argv) != 4:
# #         print("使用方式: python script.py <csv_file> <transcript_name> <transcript_length>")
# #     else:
# #         csv_file = sys.argv[1]
# #         transcript_name = sys.argv[2]
# #         transcript_length = int(sys.argv[3])
# #         main(csv_file, transcript_name, transcript_length)

# ##############
# # v2 overlap #
# ##############

# import sys
# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np

# def main(csv_file, transcript_name, transcript_length):
#     # 讀取 CSV
#     df = pd.read_csv(csv_file, comment='#')

#     # 過濾指定 transcript
#     df = df[df['ref_id'] == transcript_name]

#     if df.empty:
#         print("找不到指定的 transcript")
#         return

#     # 初始化 coverage array
#     coverage = np.zeros(transcript_length, dtype=float)

#     # 疊加每一筆 read_count 到對應區間
#     for _, row in df.iterrows():
#         start = int(row['init_pos'])
#         end = int(row['end_pos'])
#         count = float(row['read_count'])

#         # 注意避免超過 transcript 長度
#         end = min(end, transcript_length)
#         coverage[start:end + 1] += count

#     ################################################
#     ## output as csv to compare with other method ##
#     ################################################
    
#     # # 匯出 coverage 結果
#     # segments = []
#     # current_start = 0
#     # current_count = coverage[0]

#     # for i in range(1, transcript_length):
#     #     if coverage[i] != current_count:
#     #         segments.append({
#     #             'start': current_start,
#     #             'end': i-1,
#     #             'read_count': current_count
#     #         })
#     #         current_start = i
#     #         current_count = coverage[i]
#     # # 最後一段
#     # segments.append({
#     #     'start': current_start,
#     #     'end': transcript_length,
#     #     'read_count': current_count
#     # })
#     # df = pd.DataFrame(segments)
#     # df['read_count'] = df['read_count'].round(3)
#     # df[df['read_count'] != 0].to_csv("output/{}_bedgraph.csv".format(csv_file.replace(".csv", '').split('/')[-1]), index=False)
#     ################################################

#     # 準備畫圖
#     x = np.arange(transcript_length)
#     y = coverage

#     plt.figure(figsize=(12, 6))
#     plt.bar(x, y, width=1, color='skyblue', edgecolor='none')

#     # 標記最大值
#     max_y = np.max(y)
#     max_x = np.argmax(y)
#     if max_y > 0:
#         plt.text(max_x, max_y, str(max_y), fontsize=10, ha='center', va='bottom', color='red')

#     # 標籤與範圍
#     plt.xlabel('Transcript Position')
#     plt.ylabel('Read Count')
#     # plt.title('Accumulated Read Count Profile for Transcript: {}'.format(transcript_name))
#     plt.title('{} mapped to PAR-6 gene({})'.format(csv_file.replace(".csv", '').split('/')[-1], transcript_name))
#     plt.xlim(0, transcript_length)
#     plt.tight_layout()

#     # 儲存與顯示
#     output_filename = "output/{}_22G.png".format(csv_file.replace(".csv", '').split('/')[-1])
#     plt.savefig(output_filename)
#     print("save graph to:", output_filename)
#     plt.show()

# if __name__ == '__main__':
#     if len(sys.argv) != 4:
#         print("使用方式: python script.py <csv_file> <transcript_name> <transcript_length>")
#     else:
#         csv_file = sys.argv[1]
#         transcript_name = sys.argv[2]
#         transcript_length = int(sys.argv[3])
#         main(csv_file, transcript_name, transcript_length)

#########################
## 2 input same height ##
#########################

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def calculate_coverage(csv_file, transcript_name, transcript_length):
    df = pd.read_csv(csv_file, comment='#')
    df = df[df['ref_id'] == transcript_name]

    if df.empty:
        print("{}: 找不到指定的 transcript".format(csv_file))
        return None

    coverage = np.zeros(transcript_length, dtype=float)

    for _, row in df.iterrows():
        start = int(row['init_pos'])
        end = min(int(row['end_pos']), transcript_length)
        count = float(row['read_count'])
        coverage[start:end + 1] += count

    return coverage

def plot_coverage(coverage, transcript_length, transcript_name, csv_file, y_max):
    x = np.arange(transcript_length)
    y = coverage

    plt.figure(figsize=(12, 6))
    plt.bar(x, y, width=1, color='skyblue', edgecolor='none')

    max_y = np.max(y)
    max_x = np.argmax(y)
    if max_y > 0:
        plt.text(max_x, max_y, str(round(max_y, 2)), fontsize=10, ha='center', va='bottom', color='red')

    plt.xlabel('Transcript Position')
    plt.ylabel('Read Count')
    plt.title('{} mapped to PAR-6 gene({})'.format(os.path.basename(csv_file).replace(".csv", ''), transcript_name))
    plt.xlim(0, transcript_length)
    plt.ylim(0, y_max)
    plt.tight_layout()

    output_filename = "output/{}_22G.png".format(os.path.basename(csv_file).replace(".csv", ''))
    plt.savefig(output_filename)
    print("save graph to:", output_filename)
    plt.show()

def main(csv_file1, csv_file2, transcript_name, transcript_length):
    cov1 = calculate_coverage(csv_file1, transcript_name, transcript_length)
    cov2 = calculate_coverage(csv_file2, transcript_name, transcript_length)

    if cov1 is None or cov2 is None:
        return

    y_max = max(np.max(cov1), np.max(cov2)) * 1.05

    plot_coverage(cov1, transcript_length, transcript_name, csv_file1, y_max)
    plot_coverage(cov2, transcript_length, transcript_name, csv_file2, y_max)

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("使用方式: python script.py <csv_file1> <csv_file2> <transcript_name> <transcript_length>")
    else:
        csv_file1 = sys.argv[1]
        csv_file2 = sys.argv[2]
        transcript_name = sys.argv[3]
        transcript_length = int(sys.argv[4])
        main(csv_file1, csv_file2, transcript_name, transcript_length)
