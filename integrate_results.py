import pandas as pd
import numpy as np
import sys
def strip_suffix(gene_id):
    """
    遺伝子IDから転写産物のサフィックス (.t1, .1 など) を取り除くヘルパー関数。
    parse_gene_trees.py の get_gene_id と似たロジック。
    """
    # Glymaの特殊ルール
    if gene_id.startswith('Glyma.') and gene_id.count('.') >= 2:
        parts = gene_id.split('.')
        if len(parts) > 2 and parts[-1].isalpha() and parts[-2].isdigit():
            return '.'.join(parts[:-2])

    # 一般的なサフィックスルール
    parts = gene_id.split('.')
    if len(parts) > 1:
        suffix = parts[-1]
        if (suffix.startswith('t') and suffix[1:].isdigit()) or suffix.isdigit():
            return '.'.join(parts[:-1])
            
    return gene_id
    
def integrate_synteny_and_phylogeny(phylo_pairs_csv, synteny_pairs_tsv, output_csv):
    """
    系統樹ベースのオーソログペアと、MCScanXのシンテニー情報を統合し、
    最終的なマッピング判定を行う。

    Args:
        phylo_pairs_csv (str): parse_gene_trees.py の出力CSVファイル。
        synteny_pairs_tsv (str): MCScanXから作成した2列のシンテニーペアファイル。
        output_csv (str): 最終結果を出力するCSVファイル名。
    """
    # 1. データの読み込み
    print("Reading phylogenetic pairs...")
    try:
        df_phylo = pd.read_csv(phylo_pairs_csv)
    except FileNotFoundError:
        print(f"Error: Phylogenetic pairs file not found at {phylo_pairs_csv}", file=sys.stderr)
        sys.exit(1)

    print("Reading syntenic pairs...")
    try:
        # タブ区切りで、ヘッダーなしと仮定
        df_synteny = pd.read_csv(synteny_pairs_tsv, sep='\t', header=None, names=['s1_gene_synteny', 's2_gene_synteny'])
    except FileNotFoundError:
        print(f"Error: Syntenic pairs file not found at {synteny_pairs_tsv}", file=sys.stderr)
        sys.exit(1)
    
    # 高速な検索のために、シンテニーペアをセットに変換
    # (gene1, gene2) のタプルのセットを作成
    syntenic_pairs_set = set(zip(df_synteny['s1_gene_synteny'], df_synteny['s2_gene_synteny']))
    
    print(f"Loaded {len(syntenic_pairs_set)} unique syntenic pairs.")

    # 2. シンテニー情報の付加
    print("Stripping suffixes from phylogenetic IDs for comparison...")
    # 系統ペアのIDからサフィックスを除去した新しい列を一時的に作成
    df_phylo['s1_gene_stripped'] = df_phylo['s1_gene'].apply(strip_suffix)
    df_phylo['s2_gene_stripped'] = df_phylo['s2_gene'].apply(strip_suffix)

    print("Checking for syntenic relationships...")
    # サフィックス除去後のIDでシンテニー関係をチェック
    df_phylo['is_syntenic'] = df_phylo.apply(
        lambda row: (row['s1_gene_stripped'], row['s2_gene_stripped']) in syntenic_pairs_set,
        axis=1
    )
    
    # 一時的に使った列を削除して、出力ファイルをクリーンにする
    df_phylo = df_phylo.drop(columns=['s1_gene_stripped', 's2_gene_stripped'])

    # 3. 最終的なマッピング判定
    # 優先順位: 1. シンテニーあり & one-to-one > 2. シンテニーあり > 3. one-to-one > 4. その他
    conditions = [
        (df_phylo['is_syntenic'] == True) & (df_phylo['relationship'] == 'one-to-one'),
        (df_phylo['is_syntenic'] == True),
        (df_phylo['relationship'] == 'one-to-one'),
    ]
    choices = [
        'High-confidence (Syntenic & 1-to-1)',
        'High-confidence (Syntenic)',
        'Medium-confidence (1-to-1, non-syntenic)',
    ]
    
    df_phylo['mapping_status'] = np.select(conditions, choices, default='Low-confidence')

    # 4. 結果の保存
    # 判定結果でソートすると見やすい
    df_phylo_sorted = df_phylo.sort_values(by=['s1_gene', 'mapping_status'])
    df_phylo_sorted.to_csv(output_csv, index=False)
    
    print(f"\nIntegration complete. Results saved to {output_csv}")
    print("\nMapping status summary:")
    print(df_phylo_sorted['mapping_status'].value_counts())


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python integrate_results.py <phylo_pairs.csv> <synteny_pairs.tsv> <output_final.csv>")
        sys.exit(1)

    phylo_file = sys.argv[1]
    synteny_file = sys.argv[2]
    output_file = sys.argv[3]
    
    integrate_synteny_and_phylogeny(phylo_file, synteny_file, output_file)
