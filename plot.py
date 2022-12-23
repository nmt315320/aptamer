import numpy as np
import matplotlib.pyplot as plt

# x轴刻度标签
x_ticks = ['a', 'b', 'c', 'd', 'e', 'f']
# 柱的宽度
barWidth = 0.25
# 第1个柱的x轴范围（每个柱子的中点）（0, 1, ..., len(x_ticks)-1）
x1 = np.arange(len(x_ticks))
# 第2个柱的x轴范围（每个柱子的中点）
x2 = [x + barWidth for x in x1]
# 第1个柱数据
y1 = [5, 3, 2, 4, 1, 6]
# 第2个柱数据
y2 = [3, 1, 6, 5, 2, 4]

# 设置画布大小
plt.figure(figsize=(10, 6))
# 画第1个柱
plt.bar(x1, y1, color='#FF0000', width=barWidth, label='label1')
# 画第2个柱
plt.bar(x2, y2, color='#00FF00', width=barWidth, label='label2')

# 给第1个柱数据点加上数值，前两个参数是坐标，第三个是数值，ha和va分别是水平和垂直位置（数据点相对数值）。
for a, b in zip(x1, y1):
    plt.text(a, b, '%d'%b, ha='center', va= 'bottom', fontsize=18)
# 给第2个柱数据点加上数值
for a, b in zip(x2, y2):
    plt.text(a, b, '%d'%b, ha='center', va= 'bottom', fontsize=18)

# 画水平横线
plt.hlines(3, 0, len(x_ticks)-1+barWidth, colors = "#000000", linestyles = "dashed")

# 添加x轴和y轴刻度标签
plt.xticks([r + barWidth/2 for r in x1], x_ticks, fontsize=18)
plt.yticks(fontsize=18)

# 添加x轴和y轴标签0
plt.xlabel(u'x_label', fontsize=18)
plt.ylabel(u'y_label', fontsize=18)

# 标题
plt.title(u'Title', fontsize=18)

# 图例
plt.legend(fontsize=18)

# 保存图片
plt.savefig('./figure.pdf', bbox_inches='tight')
# 显示图片
plt.show()