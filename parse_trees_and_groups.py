import pandas as pd
import os
import sys
from ete3 import Tree

def find_best_ortholog_pairs_from_tree(gene_tree_path, species1_prefix, species2_prefix):
    """
    単一の遺伝子系統樹ファイルを解析し、2種間のオーソログペア候補を特定する。

    Args:
        gene_tree_path (str): OrthoFinderが出力した遺伝子系統樹ファイルのパス。
        species1_prefix (str): 比較元となる種の遺伝子IDの接頭辞 (例: 'Fvesca_v4_')。
        species2_prefix (str): 比較対象となる種の遺伝子IDの接頭辞 (例: 'Fvesca_v6_')。

    Returns:
        list: 見つかったペア情報の辞書のリスト。
              例: [{'s1_gene': 'Fvesca_v4_g1', 's2_gene': 'Fvesca_v6_gX', 
                    'relationship': 'one-to-one', 'distance': 0.123}, ...]
    """
    try:
        tree = Tree(gene_tree_path, format=1)
    except Exception as e:
        print(f"Warning: Could not parse tree {gene_tree_path}. Error: {e}", file=sys.stderr)
        return []

    # 葉（遺伝子）を種ごとに分類
    species1_leaves = [leaf for leaf in tree.iter_leaves() if leaf.name.startswith(species1_prefix)]
    species2_leaves = [leaf for leaf in tree.iter_leaves() if leaf.name.startswith(species2_prefix)]

    if not species1_leaves or not species2_leaves:
        return []

    results = []

    # species1 の各遺伝子について、最も近縁な species2 の遺伝子を探す
    for s1_leaf in species1_leaves:
        best_s2_leaf = None
        min_distance = float('inf')

        # 総当たりで最も近いペアを探す
        for s2_leaf in species2_leaves:
            # 2つの葉の間の距離を計算
            distance = s1_leaf.get_distance(s2_leaf)
            if distance < min_distance:
                min_distance = distance
                best_s2_leaf = s2_leaf
        
        if not best_s2_leaf:
            continue

        # 最も近かったペアの関係性を分類
        # 1. 共通祖先ノードを取得
        ancestor = s1_leaf.get_common_ancestor(best_s2_leaf)
        
        # 2. 共通祖先の子孫の数をチェック
        descendant_leaves = ancestor.get_leaves()
        
        s1_in_clade = [leaf for leaf in descendant_leaves if leaf.name.startswith(species1_prefix)]
        s2_in_clade = [leaf for leaf in descendant_leaves if leaf.name.startswith(species2_prefix)]
        
        relationship = 'ambiguous' # デフォルト
        if len(s1_in_clade) == 1 and len(s2_in_clade) == 1:
            relationship = 'one-to-one'
        elif len(s1_in_clade) == 1 and len(s2_in_clade) > 1:
            relationship = 'one-to-many'
        elif len(s1_in_clade) > 1 and len(s2_in_clade) == 1:
            relationship = 'many-to-one'
        elif len(s1_in_clade) > 1 and len(s2_in_clade) > 1:
            relationship = 'many-to-many'

        results.append({
            's1_gene': s1_leaf.name,
            's2_gene': best_s2_leaf.name,
            'relationship': relationship,
            'distance': min_distance,
            'orthogroup': os.path.basename(gene_tree_path).replace('_tree.txt', '')
        })
        
    return results



def find_pairs_from_group_counts(orthogroup_id, orthogroup_row, s1_col, s2_col):
    """
    【新規関数】系統樹がないオーソグループについて、遺伝子数で関係を分類する。
    """
    s1_genes_str = orthogroup_row[s1_col]
    s2_genes_str = orthogroup_row[s2_col]
    
    if pd.isna(s1_genes_str) or pd.isna(s2_genes_str):
        return []
        
    s1_genes = s1_genes_str.split(', ')
    s2_genes = s2_genes_str.split(', ')
    s1_count = len(s1_genes)
    s2_count = len(s2_genes)
    
    relationship = 'count_ambiguous'
    if s1_count == 1 and s2_count == 1: relationship = 'count_one-to-one'
    elif s1_count == 1 and s2_count > 1: relationship = 'count_one-to-many'
    elif s1_count > 1 and s2_count == 1: relationship = 'count_many-to-one'
    elif s1_count > 1 and s2_count > 1: relationship = 'count_many-to-many'
        
    return [{
        's1_gene': s1_genes[0], 's2_gene': s2_genes[0],
        'relationship': relationship, 'distance': -1.0,
        # 'orthogroup': orthogroup_row['Orthogroup'] # この行がエラーの原因
        'orthogroup': orthogroup_id  # ループ変数から受け取ったIDを使用
    }]

def main(orthofinder_results_dir, species1_id, species2_id, species1_prefix, species2_prefix, output_file):
    gene_trees_dir = os.path.join(orthofinder_results_dir, 'Gene_Trees')
    orthogroups_file = os.path.join(orthofinder_results_dir, 'Orthogroups', 'Orthogroups.tsv')

    if not os.path.isdir(gene_trees_dir) or not os.path.isfile(orthogroups_file):
        print(f"Error: Required files not found in {orthofinder_results_dir}", file=sys.stderr)
        sys.exit(1)

    # 1. 存在する系統樹のリストを作成
    existing_trees = {f.replace('_tree.txt', '') for f in os.listdir(gene_trees_dir) if f.endswith('_tree.txt')}

    # 2. Orthogroups.tsv を読み込む
    df_groups = pd.read_csv(orthogroups_file, sep='\t', index_col=0)
    
    # species IDから正しい列名を取得 (ファイル名が列名になる)
    # 例: Fvesca_v4 -> Fvesca_v4.fa.pep などの列名を探す
    s1_col_name = next((col for col in df_groups.columns if species1_id in col), None)
    s2_col_name = next((col for col in df_groups.columns if species2_id in col), None)

    if not s1_col_name or not s2_col_name:
        print(f"Error: Could not find columns for '{species1_id}' or '{species2_id}' in Orthogroups.tsv", file=sys.stderr)
        print("Available columns:", df_groups.columns.tolist())
        sys.exit(1)
    
    print(f"Found columns: {s1_col_name} for species1, {s2_col_name} for species2")

    all_pairs = []
    
    # 3. 全オーソグループをループ処理
    print(f"Processing {len(df_groups)} orthogroups...")
    for orthogroup, row in df_groups.iterrows():
        # 3a. 系統樹が存在する場合
        if orthogroup in existing_trees:
            tree_path = os.path.join(gene_trees_dir, f"{orthogroup}_tree.txt")
            pairs = find_best_ortholog_pairs_from_tree(tree_path, species1_prefix, species2_prefix)
            all_pairs.extend(pairs)
        # 3b. 系統樹が存在しない場合
        else:
            pairs = find_pairs_from_group_counts(orthogroup, row, s1_col_name, s2_col_name)
            all_pairs.extend(pairs)

    print(f"Total pairs/representatives found: {len(all_pairs)}")

    # 4. 結果をCSVファイルに書き出し
    with open(output_file, 'w') as f_out:
        f_out.write("s1_gene,s2_gene,relationship,distance,orthogroup\n")
        for pair in all_pairs:
            f_out.write(f"{pair['s1_gene']},{pair['s2_gene']},{pair['relationship']},{pair['distance']},{pair['orthogroup']}\n")

    print(f"Results successfully written to {output_file}")

if __name__ == '__main__':
    if len(sys.argv) != 7:
        print("Usage: python parse_trees_and_groups.py <orthofinder_results_dir> <species1_id> <species2_id> <species1_prefix> <species2_prefix> <output_file.csv>")
        print("\nExample:")
        print("python parse_trees_and_groups.py orthofinder_results/ 'Fvesca_v4' 'Fvesca_v6' 'v4' 'v6' v4_v6_all_orthologs.csv")
        sys.exit(1)

    orthofinder_dir = sys.argv[1]
    s1_id = sys.argv[2]
    s2_id = sys.argv[3]
    s1_prx = sys.argv[4]
    s2_prx = sys.argv[5]
    output_csv = sys.argv[6]

    main(orthofinder_dir, s1_id, s2_id, s1_prx, s2_prx, output_csv)
