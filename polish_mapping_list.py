import pandas as pd
import sys

def polish_mapping_list(input_csv, output_csv):
    """
    最終マッピングリストを磨き上げる。
    信頼性の低いペアのうち、s2遺伝子が既に信頼性の高いペアで
    使用されているものを削除する。
    """
    # 1. データの読み込み
    print(f"Reading final mapping list from {input_csv}...")
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_csv}", file=sys.stderr)
        sys.exit(1)
        
    print(f"Total pairs before polishing: {len(df)}")

    # 2. 信頼度に基づいてペアをグループ分け
    # "High" または "Medium" を含むものを Confident とする
    is_confident = df['mapping_status'].str.contains('High-confidence|Medium-confidence', case=False, na=False)
    df_confident = df[is_confident]
    df_low_confidence = df[~is_confident]

    print(f"Found {len(df_confident)} confident pairs.")
    print(f"Found {len(df_low_confidence)} low-confidence pairs to check.")

    # 3. Confidentなペアで使われているs2遺伝子のセットを作成
    confident_s2_genes = set(df_confident['s2_gene'])
    print(f"Found {len(confident_s2_genes)} unique s2_genes in confident pairs.")

    # 4. Low-confidenceペアをフィルタリング
    # s2遺伝子がconfident_s2_genesセットに含まれていないペアだけを残す
    df_low_confidence_filtered = df_low_confidence[~df_low_confidence['s2_gene'].isin(confident_s2_genes)]
    
    num_removed = len(df_low_confidence) - len(df_low_confidence_filtered)
    print(f"Removed {num_removed} low-confidence pairs whose s2_gene was already confidently mapped.")
    print(f"Kept {len(df_low_confidence_filtered)} low-confidence pairs that have no confident partner.")

    # 5. フィルタリング後のリストを結合
    df_final_polished = pd.concat([df_confident, df_low_confidence_filtered], ignore_index=True)
    
    # 見やすいようにソート
    df_final_polished = df_final_polished.sort_values(by='s1_gene')
    
    # 6. 結果を保存
    df_final_polished.to_csv(output_csv, index=False)

    print(f"\nPolishing complete. Total pairs in final list: {len(df_final_polished)}")
    print(f"Polished mapping list saved to {output_csv}")
    print("\nFinal mapping status summary:")
    print(df_final_polished['mapping_status'].value_counts())


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python polish_mapping_list.py <input_final_mapping.csv> <output_polished_mapping.csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    polish_mapping_list(input_file, output_file)
