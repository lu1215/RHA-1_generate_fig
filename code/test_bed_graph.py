import pandas as pd
import numpy as np
import sys
import os

def process_single_transcript(input_csv, output_csv, transcript_name='T26E3.3a.1'):
    # 讀取 CSV 並過濾指定 transcript
    df = pd.read_csv(input_csv, comment='#')
    df = df[df['ref_id'] == transcript_name].sort_values(by='init_pos').reset_index(drop=True)

    if df.empty:
        print("找不到 transcript: {}".format(transcript_name))
        return

    # 初始化 coverage array
    bottom = df['init_pos'].min()
    top = df['end_pos'].max()
    arr = np.zeros(top - bottom + 1)

    # 累加 read_count 到每個位置
    for i in range(len(df)):
        start = df.loc[i, 'init_pos']
        end = df.loc[i, 'end_pos']
        count = df.loc[i, 'read_count']
        arr[start - bottom : end - bottom + 1] += count

    # 壓縮連續相同 read_count 的區段
    init_pos = []
    end_pos = []
    read_count = []

    for i in range(len(arr)):
        if i == 0:
            init_pos.append(bottom)
        if i == len(arr) - 1:
            end_pos.append(bottom + i)
            read_count.append(arr[i])
            break
        if arr[i] != arr[i+1]:
            end_pos.append(bottom + i)
            init_pos.append(bottom + i + 1)
            read_count.append(arr[i])

    # 組成 DataFrame 並重新命名欄位順序
    result_df = pd.DataFrame({
        'start': init_pos,
        'end': end_pos,
        'read_count': read_count
    })

    # 過濾掉 read_count == 0 的區段
    result_df = result_df[result_df['read_count'] != 0]

    # 調整欄位順序為 end, read_count, start（你指定的格式）
    result_df = result_df[['end', 'read_count', 'start']]
    result_df['read_count'] = result_df['read_count'].round(3)

    # 確保資料照 start 排序
    result_df = result_df.sort_values(by='start').reset_index(drop=True)

    # 初始化合併清單
    merged_rows = []

    # 初始化第一段
    current_start = result_df.loc[0, 'start']
    current_end = result_df.loc[0, 'end']
    current_count = result_df.loc[0, 'read_count']

    for i in range(1, len(result_df)):
        row = result_df.loc[i]
        if row['read_count'] == current_count and row['start'] == current_end + 1:
            # 可以合併 → 延長 end
            current_end = row['end']
        else:
            # 儲存目前段落
            merged_rows.append({
                'start': int(current_start),
                'end': int(current_end),
                'read_count': current_count
            })
            # 開始新段落
            current_start = row['start']
            current_end = row['end']
            current_count = row['read_count']

    # 別忘了最後一段
    merged_rows.append({
        'start': int(current_start),
        'end': int(current_end),
        'read_count': current_count
    })

    # 轉回 DataFrame，並照原順序：end, read_count, start
    result_df = pd.DataFrame(merged_rows)[['end', 'read_count', 'start']]

    # 建立輸出資料夾（相容 Python 3.5）
    out_dir = os.path.dirname(output_csv)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # 輸出結果
    result_df.to_csv(output_csv, index=False)
    print("結果已儲存為：{}".format(output_csv))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("使用方式: python test_bed_graph.py <input_csv> <output_csv>")
    else:
        input_csv = sys.argv[1]
        output_csv = sys.argv[2]
        process_single_transcript(input_csv, output_csv)
