# 提交要求

每个提交（包括仅改文档）都必须在 `hzcu-hpct` 通过 `make bench-threads` 的全线程矩阵：1/2/4/8/16/32/64 线程均绑定在 NUMA 1，64 线程使用 32 个物理核的 SMT。每档 5 次，逐次核对官方输出。不得用 `--no-verify` 跳过钩子，测试失败不得提交。

新检出须先设置 `git config --local core.hooksPath .githooks`。钩子测试暂存树；不触发钩子的 Git 操作也须逐提交复测。资源检查、日志位置与复现步骤见 `reports/README.md`。
