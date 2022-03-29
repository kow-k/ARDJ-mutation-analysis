# ARDJ-mutation-analysis

R scripts and data used in my D3-3 paper presented at NLP28/2022

2022年3月に実施された言語処理学会28回年次大会 (NLP28/2022) のD3-3の発表 "How nearly random mutations to sentences affect their acceptabilities: A preliminary analysis" (Kow Kuroda, Hikaru Yokono, Keiga Abe, Tomoyuki Tsuchiya, Yoshihiko Asao, Yuichiro Kobayashi, Toshiyuki Kanamaru and Takumi Tagawa) で使った解析処理とデータの公開

# R scripts

1, 2, 3, 4 の順に実行する:

1. [cluster-originals.R](cluster-originals.R) [前処理と 36原文の DBSCAN クラスタリング]

2. [analyze-mutations-by-cluster.R](analyze-mutations-by-clusters.R) [DBSCAN クラスター1 の分析]

3. [cluster-originals-by-KL-divergence.R](cluster-originals-by-KL-divergence.R) [36原文の KL-divergence 分析]

4. [analyze-by-KL-divergence.R](analyze-by-KL-divergence.R) [DBSCAN クラスター1 の KL-divergence 分析]



# data

上の R scripts が処理する [応答データ](data-s2u-sd-filtered1.csv)


# paper

発表した [論文](https://www.dropbox.com/s/ziesrj8o8yrz49b/ARDJ-mutation-analysis-nlp28.pdf?dl=0)

# slide

発表に使った[スライド (付録つき)](https://www.dropbox.com/s/atupwiesjqbzlss/ARDJ-mutation-analysis-nlp28-slides.pdf?dl=0)


