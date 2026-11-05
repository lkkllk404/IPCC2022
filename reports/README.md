# IPCC2022 优化记录

完成 **5 个性能提交**。统一使用官方 `uniformvector-2dim-5h.txt`（dim=2,n=500,k=2），gcc12.3.0 `-O3 -march=native`，Xeon Gold 6430 上单逻辑核 CPU100，三次 wall time 中位数。每轮详细记录涵盖上次提交以来的全部候选、测试、资源估算和结果；没有安装依赖或修改服务器环境。

|阶段|变更|中位数(s)|较上一版本|累计加速|报告|
|---|---|---:|---:|---:|---|
|基线 c3a01b9|仅修复重复 fclose 的 UB|61.70|—|1.00x|[01](01.md)|
|1 / 8c47585|点对对称性|28.04|2.20x|2.20x|[01](01.md)|
|2 / 73e7062|距离缓存、k=2 连续访存|12.19|2.30x|5.06x|[02](02.md)|
|3 / 185d9ca|Chebyshev→两个一维投影排序|9.79|1.25x|6.30x|[03](03.md)|
|4 / 678a4d5|类型专用比较排序|6.03|1.62x|10.23x|[04](04.md)|
|5 / 本提交|有界字节基数排序|2.51|2.40x|24.58x|[05](05.md)|

## 复现

每次开始测试前先检查 `uptime`、`free -m` 和其他任务占用。默认只用一个逻辑CPU，官方输入最终版本RSS约4.1 MiB；不要直接改成多核/大规模输入而不先评估资源。

```sh
cd ~/Sandboxes/HPC/IPCC2022
uptime; free -m
taskset -c 100 make -j1
make check        # 使用已有 /opt/anaconda3/bin/python；无需第三方包
make bench        # 串行3次、CPU100、timeout180秒，严格cmp官方输出

# 可选择重跑初始修复版，约62秒/次；不要与其他测试并发
# taskset -c 100 make reports/work/baseline
# bash reports/benchmark.sh reports/work/baseline baseline 3
```

`BENCH_CPU=<cpu>` 可更改绑核；`make PYTHON=<已有conda解释器>` 可更改 Python 路径。`-march=native` 产物仅面向构建机器，换机器需重新编译。`reports/work/` 忽略临时产物；每轮的 `.time`、验证日志和摘要已经纳入 Git。

## 保证范围

- 所有轮次官方输出与 `refer-2dim-5h.txt` 逐字节一致。
- 末轮通过10个端到端场景、945个直接分数测试、112个排序场景、11个错误处理场景及ASan/UBSan。
- 对浮点等分，组合之间可能换位；不是任意输入逐字节复刻原版的保证。差分测试会检查每个名次独立重算的数学分数，官方输入则必须严格cmp。
- k=2 的单组合复杂度为 O(n log n)，在基数排序适用范围内近似 O(n)；k!=2 保留通用点对算法。枚举组合的数量仍是 C(n,k)，不承诺解决任意大规模输入。
- 通用缓存内存 O(n²+nk+Mk)，大 n 运行前必须重新估算，不应照搬本次输入的低内存数字。
