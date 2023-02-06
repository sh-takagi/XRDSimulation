# XRDSimulation

### できること
XRDの実験データとその構造のcifファイルから結晶の長さ、角度を求める。

### 方法
インプットのファイルを用意して同じファイル内で`python XRDSimulation.py [trials]`を実行  
`[trials]`は計算回数を整数型で入力

### インプット
- XRDの実験データ(filenameをtarget.csvにする)
- 構造のcifファイル(filenameをstructure.cifにする)

### アウトプット
Best Score  
{   
  a: a軸の長さ  
  b: b軸の長さ  
  c: c軸の長さ  
  alpha: 角度alpha  
  beta:  角度beta  
  gamma: 角度gamma  
}
