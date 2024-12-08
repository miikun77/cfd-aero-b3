# cfd-aero-b3

## 動かし方

```bash
python3 main.py
```

もしくは

```bash
python main.py
```

## 変更すべき箇所

- `t_mesh.f` の `CDEFALULT`
- `mesh.d` の 各種数値 (`main.py` で指定)
- `airfoil.dat` をどのようなファイルにするか（初期）

## 生成ファイル

`./work` にあります。計算に使用しているものは全てここ。
計算結果は`result.csv`にあります。
