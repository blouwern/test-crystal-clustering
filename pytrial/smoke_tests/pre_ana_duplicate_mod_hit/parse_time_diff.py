import re
import matplotlib.pyplot as plt


def parse_file(filename):
    consec_diffs = []
    first_diffs = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            # 只处理包含目标关键字的行
            if "Consecutive diffs" in line or "Diffs from 1st hit" in line:
                # 提取方括号内容
                match = re.search(r"\[([^\]]*)\]", line)
                if match:
                    nums_str = match.group(1)  # 括号内的字符串
                    # 提取数字（支持整数、小数、正负号）
                    nums = re.findall(r"[-+]?\d*\.?\d+", nums_str)
                    diffs = [float(x) for x in nums if x]  # 过滤空字符串
                    if "Consecutive" in line:
                        consec_diffs.extend(diffs)
                    else:
                        first_diffs.extend(diffs)
                else:
                    # 如果实际文件没有方括号（备选），可以用下面这一行提取整行数字，
                    # 但可能误提取 "1st" 中的 1，所以不推荐。
                    # 此处加个警告，让你知道格式不符
                    print(f"Warning: No brackets found in line: {line}")
    return consec_diffs, first_diffs


consec, first = parse_file("data_duplicate_modid_total.log")
print(consec)
print(first)

# 绘制直方图
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.hist(consec, bins=10, range=(0, 120), color="blue", alpha=0.7)
ax1.set_title("Consecutive diffs")
ax1.set_xlabel("Delta t (ns)")
ax1.set_ylabel("Entries")

ax2.hist(first, bins=10, range=(0, 120), color="red", alpha=0.7)
ax2.set_title("Diffs from 1st hit")
ax2.set_xlabel("Delta t (ns)")
ax2.set_ylabel("Entries")

plt.tight_layout()
plt.savefig("time_diffs.png")
plt.show()
