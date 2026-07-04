# 微分方程数值解法

李荣华 李永海 武海军 编著

# 内容提要

本书是作者结合多年来的教学经验，为适应新时代教学和发展而编写的。全书共十一章，第1章介绍常微分方程初值问题的数值解法；第2、3、4章分别探讨椭圆型、抛物型和双曲型方程的有限差分法；第5到9章深入讨论边值问题的变分形式与Ritz-Galerkin法、有限元法及其多种变体，包括有限体积元法、间断Galerkin法和弱有限元法；第10、11章介绍有限元多重网格法和自适应算法。本书专为信息与计算科学专业本科生设计，同时也适用于应用数学、力学及工程科学专业的教学，并可为从事科学技术及工程计算的专业人员提供参考。

# 总序

自数学出现以来，世界上不同国家、地区的人们在生产实践中、在思考探索中以不同的节奏推动着数学的不断突破和飞跃，并使之成为一门系统的学科。尤其是进入21世纪之后，数学发展的速度、规模、抽象程度及其应用的广泛和深入都远远超过了以往任何时期。数学的发展不仅是在理论知识方面的增加和扩大，更是思维能力的转变和升级，数学深刻地改变了人类认识和改造世界的方式。对于新时代的数学研究和教育工作者而言，有责任将这些知识和能力的发展与革新及时体现到课程和教材改革等工作当中。

数学“101计划”核心教材是我国高等教育领域数学教材的大型编写工程。作为教育部基础学科系列“101计划”的一部分，数学“101计划”旨在通过深化课程、教材改革，探索培养具有国际视野的数学拔尖创新人才，教材的编写是其中一项重要工作。教材是学生理解和掌握数学的主要载体，教材质量的高低对数学教育的变革与发展意义重大。优秀的数学教材可以为青年学生打下坚实的数学基础，培养他们的逻辑思维能力和解决问题的能力，激发他们进一步探索数学的兴趣和热情。为此，数学“101计划”工作组统筹协调来自国内16所一流高校的师资力量，全面梳理知识点，强化协同创新，陆续编写完成符合数学学科“教与学”特点，体现学术前沿，具备中国特色的高质量核心教材。此次核心教材的编写者均为具有丰富教学成果和教材编写经验的数学家，他们当中很多人不仅有国际视野，还在各自的研究领域作出杰出的工作成果。在教材的内容方面，几乎是包括了分析学、代数学、几何学、微分方程、概率论、现代分析、数论基础、代数几何基础、拓扑学、微分几何、应用数学基础、统计学基础等现代数学的全部分支方向。考虑到不同层次的学生需要，编写组对个别教材设置了不同难度的版本。同时，还及时结合现代科技的最新动向，特别组织编写《人工智能的数学基础》等相关教材。

数学“101计划”核心教材得以顺利完成离不开所有参与教材编写和审订的专家、学者及编辑人员的辛勤付出，在此深表感谢。希望读者们能通过数学“101计划”核心教材更好地构建扎实的数学知识基础，锻炼数学思维能力，深化对数学的

理解，进一步生发出自主学习探究的能力。期盼广大青年学生受益于这套核心教材，有更多的拔尖创新人才脱颖而出！

田刚

数学“101计划”工作组组长

中国科学院院士

北京大学讲席教授

# 前言

微分方程数值解法在数值分析领域占据着举足轻重的地位，它以逼近论和数值代数等学科为基础，并反过来推动这些学科的发展。微分方程数值解法广泛应用于物理、工程、经济、生物学等领域，成为多学科交叉的纽带。自20世纪40年代以来，微分方程数值解法已发展成一个庞大的计算技术学科，并成为信息与计算科学专业的基础课程之一。

为了培养拔尖创新人才，推动教育教学系统改革，笔者根据数学领域“101计划”核心教材建设的整体思路以及“大师引领、科教融汇、融通中外”的建设模式，在充分借鉴国内外优秀教材建设经验的基础上，编写了这本《微分方程数值解法》。本教材继承了李荣华与刘播编写的《微分方程数值解法(第四版)》一书中的常微分方程数值解法、椭圆型方程有限差分法、抛物型方程有限差分法、双曲型方程有限差分法和偏微分方程Ritz-Galerkin法等内容。同时，改写了有限元法相关章节。另外，为了反映微分方程数值解法的学科进展，本教材还增加了有限体积元法、间断Galerkin法、弱有限元法、多重网格法和自适应有限元法等章节。

在选材上，本书遵循少而精和可接受性强的原则，力求选取基本且对本学科发展有重要影响的内容。由于我国开设信息与计算科学专业的院校众多，各院校的情况和要求差异较大，因此教师在讲授时可根据具体情况适当删减部分内容。除文中打星号章节外，其他部分也可进一步精简，但要确保所传授知识的系统性。

虽然编者在编写过程中付出了大量努力，但书中仍可能存在一些缺点甚至错误，恳请广大师生和读者指正。在成书过程中，浙江大学包刚院士和吉林大学张然教授在整体框架和内容选择方面提出了许多宝贵建议；吉林大学吕俊良教授在教材编写中付出了大量精力；吉林大学翟起龙教授为弱有限元法一章的编写提供了诸多帮助；高等教育出版社编辑们为本书的出版付出了大量辛勤劳动。在此，向他们表示衷心的感谢！

编者

2024年8月

![](images/ad1a7dd096343baacca1c2a25d59554b69072179e40ab80e0046952d3318a8a8.jpg)

# 目录

# 第1章 常微分方程初值问题的数值解法 1

1.1 引论 2

1.1.1 一阶常微分方程初值问题 2  
1.1.2 Euler法 2  
1.1.3 线性差分方程 6  
1.1.4 Gronwall 不等式 9  
1.1.5 习题 10

1.2 线性多步法 11

1.2.1 数值积分法 11  
1.2.2 待定系数法 17  
1.2.3 预估-校正算法 19  
1.2.4 多步法的计算问题 21  
1.2.5 习题 21

1.3 相容性、稳定性和误差估计 22

1.3.1 局部截断误差和相容性 22  
1.3.2 稳定性 23  
1.3.3 收敛性和误差估计 27  
1.3.4 习题 28

1.4 单步法和 Runge-Kutta 法 29

1.4.1 Taylor展开法 29  
1.4.2 单步法的稳定性和收敛性 30  
1.4.3 Runge-Kutta法 32  
1.4.4 习题 36

*1.5 绝对稳定性和绝对稳定域 37

1.5.1 绝对稳定性 37  
1.5.2 绝对稳定域 39  
1.5.3 应用例子 40

1.5.4 习题 43

*1.6 一阶方程组和刚性问题 43

1.6.1 对一阶方程组的推广 43  
1.6.2 刚性问题 45  
1.6.3 A稳定性 46  
1.6.4 数值例子 48

*1.7 外推法 50

1.7.1 多项式外推 50  
1.7.2 对初值问题的应用 51  
1.7.3 用外推法估计误差 52   
1.7.4 习题 53

# 第2章 椭圆型方程的有限差分法 55

2.1 差分逼近的基本概念 56  
2.2 一维差分格式 60

2.2.1 直接差分化 60  
2.2.2 有限体积法（积分插值法） 63  
2.2.3 边值条件的处理 65  
2.2.4 习题 66

2.3 矩形网的差分格式 66

2.3.1 五点差分格式 67  
2.3.2 边值条件的处理 70  
2.3.3 习题 72

2.4 三角网的差分格式 73

2.4.1 习题 77

2.5 极值定理和敛速估计 77

2.5.1 差分方程 77  
2.5.2 极值定理 79  
2.5.3 五点差分格式的敛速估计 81  
2.5.4 习题 82

# 第3章 抛物型方程的有限差分法 83

3.1 最简差分格式 84  
3.1.1 习题 89

3.2 稳定性与收敛性 90

3.2.1 稳定性概念 90   
3.2.2 判别稳定性的直接估计法（矩阵法） 92

3.2.3 收敛性与敛速估计 95  
3.2.4 习题 97

3.3 Fourier 方法 97

3.3.1 习题 103

*3.4 判别差分格式稳定性的代数准则 103  
3.4.1 习题 8.8.3 108

第4章 双曲型方程的有限差分法 109

4.1 波动方程的差分逼近 110

4.1.1 波动方程及其特征 110  
4.1.2 显格式 111  
4.1.3 稳定性分析 113  
4.1.4 隐格式 116  
4.1.5 数值例子 117  
4.1.6 习题 118

4.2 一阶线性双曲方程组 118

4.2.1 双曲型方程组及其特征 118  
4.2.2 Cauchy 问题、依存域、影响域和决定域 122  
4.2.3 初边值问题 123  
4.2.4 习题 124

4.3 初值问题的差分逼近 125

4.3.1 迎风格式 125  
4.3.2 积分守恒差分格式 128  
4.3.3 粘性差分格式 130  
4.3.4 其他差分格式 132  
4.3.5 习题 133

*4.4 初边值问题和对流占优扩散方程 134

4.4.1 初边值问题 134  
4.4.2 对流占优扩散方程 136  
4.4.3 数值例子 138  
4.4.4 习题 139

第5章 边值问题的变分形式与Ritz-Galerkin法141

5.1 二次函数的极值 142  
5.1.1 习题 144  
5.2 Sobolev空间初步 144   
5.2.1 Sobolev空间的定义 144

5.2.2 Sobolev空间的性质 146   
5.2.3 习题 150

# 5.3 两点边值问题 8 150

5.3.1 极小位能原理 150  
5.3.2 虚功原理 154  
5.3.3 习题 155

# 5.4 二阶椭圆边值问题 156

5.4.1 极小位能原理 156  
5.4.2 自然边值条件 159  
5.4.3 虚功原理 160  
5.4.4 习题 161

# 5.5 Ritz-Galerkin 方法 162

5.5.1 习题 167

# *5.6 谱方法 167

5.6.1 三角函数逼近 167  
5.6.2 Fourier 谱方法 170  
5.6.3 拟谱方法（配置法） 174

# 第6章有限元法 177

# 6.1 一维例子 178

6.1.1 两点边值问题及其变分公式 178  
6.1.2 有限元方法 179  
6.1.3 有限元方程组 179  
6.1.4 先验误差估计 181

# 6.2 有限元空间的构造 183

6.2.1 有限元及有限元空间 183  
6.2.2 一维高次元 186  
6.2.3 解二维问题的矩形元与四边形元 192  
6.2.4 三角形元 197  
6.2.5 三维有限元 201  
6.2.6 习题 203

# 6.3 二阶椭圆型方程的有限元法 204

6.3.1 有限元离散 204  
6.3.2 有限元方程组的形成 205  
6.3.3 习题 209

$^{*}6.4$ 有限元法的收敛性理论 209

6.4.1 插值理论 210  
6.4.2 误差估计 216   
6.4.3 习题 217

6.5 初边值问题的有限元法 218

6.5.1 热传导方程 218  
6.5.2 波动方程 220

# 第7章 有限体积元法 223

7.1 三角形网格上有限体积元法 224

7.1.1 试探函数空间和检验函数空间 224  
7.1.2 线性元有限体积法 226  
7.1.3 稳定性分析 226  
7.1.4 误差估计 230

7.2 四边形网格上的有限体积元法 233

7.2.1 四边形网格剖分及对偶剖分 233  
7.2.2 试探函数空间和检验函数空间 235  
7.2.3 等参双线性有限体积元法 236  
7.2.4 收敛性 237  
7.2.5 数值算例 238  
7.2.6 习题 239

# 第8章 间断 Galerkin法 241

8.1 内罚间断 Galerkin 法 242

8.1.1 离散格式 242  
8.1.2 对称内罚间断 Galerkin (SIPG) 法的误差分析 244  
8.1.3 非对称内罚间断 Galerkin (NIPG) 法的误差分析 247

8.2 局部间断 Galerkin (LDG) 法 249

8.2.1 离散格式 249  
8.2.2 原始变量形式 250  
8.2.3 误差估计 252

8.3 杂交间断 Galerkin (HDG) 法 256

8.3.1 离散格式 256  
8.3.2 变分形式I 258   
8.3.3 误差估计 259   
8.3.4 后处理 262

8.3.5 变分形式II 264   
8.3.6 习题 265

# 第9章 弱有限元法 267

9.1 弱微分算子 268

9.1.1 广义弱微分算子 268  
9.1.2 离散弱微分算子 269

9.2 弱有限元数值格式 269

9.2.1 适定性 原理 光学制 食 270  
9.2.2 $L^2$ 投影的误差估计 271  
9.2.3 $H^{1}$ 误差估计 272   
9.2.4 $L^2$ 误差估计 275

9.3 无稳定子弱有限元数值格式 277

9.3.1 适定性 278  
9.3.2 $H^{1}$ 误差估计 279   
9.3.3 $L^2$ 误差估计 282

# 第10章 有限元多重网格法 285

10.1 模型问题 286  
10.2 经典迭代法 287

10.2.1 矩阵形式和算子形式 287  
10.2.2 磨光性质 289

10.3 多重网格 $V$ 循环算法 290  
10.4 完全多重网格法和工作量估计 293  
10.5 多重网格 $V$ 循环算法的矩阵形式 295

10.5.1 习题 296

# 第11章 自适应有限元法 297

11.1 一个带奇性的例子 298  
11.2 后验误差分析 300

11.2.1 Scott-Zhang插值算子 300  
11.2.2 后验误差估计 302

11.3 自适应算法 305  
11.4 收敛性分析 306

11.4.1 习题 309

# 参考文献 310

# 数学家简介 312

# 常微分方程初值问题的数值解法

# 1.1 引论

# 1.1.1 一阶常微分方程初值问题

设 $f(t,u)$ 在区域 $G:0\leqslant t\leqslant T$ ， $|u| < \infty$ 上连续，求 $u = u(t)$ 满足

$$
\frac {\mathrm {d} u}{\mathrm {d} t} = f (t, u), \quad 0 <   t \leqslant T, \tag {1.1a}
$$

$$
u (0) = u _ {0}, \tag {1.1b}
$$

其中 $u_{0}$ 是给定的初值, 这就是一阶常微分方程的初值问题. 为使问题 (1.1a) — (1.1b) 的解存在、唯一且连续依赖初值 $u_{0}$ , 即初值问题 (1.1a) — (1.1b) 适定, 还必须对右端 $f(t,u)$ 加适当限制, 通常要求 $f$ 关于 $u$ 满足 Lipschitz 条件: 存在常数 $L$ , 使

$$
\left| f \left(t, u _ {1}\right) - f \left(t, u _ {2}\right) \right| \leqslant L \left| u _ {1} - u _ {2} \right| \tag {1.2}
$$

对所有 $t\in [0,T]$ 和 $u_{1},u_{2}\in (-\infty , + \infty)$ 成立(参见[2]).本章总假定 $f$ 满足上述条件.

虽然初值问题 (1.1a) — (1.1b) 对很大一类右端函数有解, 但求出所需的解绝非易事. 事实上, 除了极特殊情形外, 人们不可能求出它的精确解, 只能用各种近似方法得到满足一定精度的近似解. 读者在常微分方程教程中已经熟悉了级数解法和 Picard 逐步逼近法, 这些方法可以给出解的近似表达式, 称为近似解析方法. 另一类近似方法只给出解在一些离散点上的近似值, 称为数值方法. 由于后一类方法应用范围更广, 特别适合用计算机计算, 所以本章只讨论初值问题的数值解法.

# 1.1.2 Euler法

最简单的数值解法是 Euler 法. 将区间 $[0, T]$ 作 $N$ 等分, 小区间的长度 $h = \frac{T}{N}$ 称为步长, 点列 $t_n = nh$ ( $n = 0, 1, \dots, N$ ) 称为节点, $t_0 = 0$ . 由已知初值 $u(t_0) = u_0$ , 可算出 $u(t)$ 在 $t = t_0$ 的导数值 $u'(t_0) = f(t_0, u(t_0)) = f(t_0, u_0)$ . 利用 Taylor 展式

$$
\begin{array}{l} u \left(t _ {1}\right) = u \left(t _ {0} + h\right) = u \left(t _ {0}\right) + h u ^ {\prime} \left(t _ {0}\right) + \frac {h ^ {2}}{2} u ^ {\prime \prime} \left(t _ {0}\right) + \frac {h ^ {3}}{6} u ^ {\prime \prime \prime} (\zeta) \\ = u _ {0} + h f \left(t _ {0}, u _ {0}\right) + R _ {0}, \tag {1.3} \\ \end{array}
$$

其中 $\zeta \in (t_0, t_1)$ . 略去二阶小量 $R_0$ , 得

$$
u _ {1} = u _ {0} + h f \left(t _ {0}, u _ {0}\right).
$$

$u_{1}$ 就是 $u(t_{1})$ 的近似值. 利用 $u_{1}$ 又可算出 $u_{2}$ , 如此下去可算出 $u$ 在所有节点上的近似

值, 一般递推公式为

$$
u _ {n + 1} = u _ {n} + h f \left(t _ {n}, u _ {n}\right), \quad n = 0, 1, \dots , N - 1,
$$

这就是Euler法

Euler法有明显的几何意义.实际上，(1.1a）的解是 $tu$ 平面上的积分曲线族，过任一点恰有一积分曲线通过.按Euler法，过初始点 $(t_0,u_0)$ 作经过此点的积分曲线的切线（斜率为 $f(t_{0},u_{0})$ )，沿切线取点 $(t_{1},u_{1})$ （ $u_{1}$ 按(1.4)计算）作为 $(t_{1},u(t_{1}))$ 的近似；然后过 $(t_{1},u_{1})$ 作一经过此点的积分曲线的切线，沿切线取点 $(t_{2},u_{2})$ （ $u_{2}$ 按(1.4)计算）作为 $(t_{2},u(t_{2}))$ 的近似.如此下去.即得一以 $(t_n,u_n)$ 为顶点的折线，这就是用Euler法得到的近似积分曲线(图1.1中的虚折线).从几何上看， $h$ 越小，此折线逼近积分曲线越好，因此也称Euler法为Euler折线法.

![](images/b5497ad610a2d1177818081a3e866f241f4136e2189d3d4d75a8381e83f218df.jpg)  
图1.1 近似积分曲线

现在用数值积分法推导 Euler 法. 将问题 (1.1a) 一 (1.1b) 写成等价的积分形式:

$$
u (t) = u _ {0} + \int_ {t _ {0}} ^ {t} f (\tau , u (\tau)) \mathrm {d} \tau \quad (t _ {0} = 0). \tag {1.5}
$$

特别地,

$$
u \left(t _ {1}\right) = u _ {0} + \int_ {t _ {0}} ^ {t _ {1}} f (\tau , u (\tau)) \mathrm {d} \tau \quad \left(t _ {0} = 0\right).
$$

用左矩形公式近似右端积分, 并用 $u_{1}$ 代替 $u(t_{1})$ 即得 $u_{1} = u_{0} + hf(t_{0},u_{0})$ , 这就是 Euler 法 (1.4). 我们也可用梯形公式近似上述积分, 仍用 $u_{1}$ 替代 $u(t_{1})$ , 得

$$
u _ {1} = u _ {0} + \frac {h}{2} [ f (t _ {0}, u _ {0}) + f (t _ {1}, u _ {1}) ].
$$

一般而言，

$$
u _ {n + 1} = u _ {n} + \frac {h}{2} [ f (t _ {n}, u _ {n}) + f (t _ {n + 1}, u _ {n + 1}) ], \quad n = 0, 1, \dots , N - 1, \tag {1.6}
$$

称之为改进的 Euler 法. 显然改进的 Euler 法比 Euler 法精度更高, 但每步计算要解非线性方程 (1.6) (关于 $u_{n+1}$ ), 这可用如下迭代公式:

$$
u _ {n + 1} ^ {[ k + 1 ]} = u _ {n} + \frac {h}{2} \left[ f \left(t _ {n}, u _ {n}\right) + f \left(t _ {n + 1}, u _ {n + 1} ^ {[ k ]}\right) \right], \quad k = 0, 1, \dots . \tag {1.7}
$$

取初值为 $u_{n + 1}^{[0]} = u_n$ ，一般只需迭代几步即可收敛

现在分析一下Euler法误差的来源.为使问题简化，我们不考虑因计算机字长限制引起的舍入误差.注意(1.3)或其一般的递推式

$$
u \left(t _ {n + 1}\right) = u \left(t _ {n}\right) + h f \left(t _ {n}, u \left(t _ {n}\right)\right) + R _ {n} \tag {1.8}
$$

是精确方程, 其中

$$
R _ {n} = \frac {h ^ {2}}{2} u ^ {\prime \prime} \left(t _ {n}\right) + \frac {h ^ {3}}{6} u ^ {\prime \prime \prime} (\zeta), \quad \zeta \in \left(t _ {n}, t _ {n + 1}\right). \tag {1.9}
$$

由(1.8)到Euler法(1.4)的唯一差别是舍去了余项 $R_{n}$ . 令

$$
L \left[ u _ {n}; h \right] = u _ {n + 1} - u _ {n} - h f \left(t _ {n}, u _ {n}\right), \tag {1.10}
$$

取 $u_{n} = u\left(t_{n}\right)$ ，则 $R_{n} = L\left[u\left(t_{n}\right);h\right] = u\left(t_{n + 1}\right) - u\left(t_{n}\right) - h u^{\prime}\left(t_{n}\right)$ .今后称 $R_{n}$ 为局部截断误差.显然Euler法的局部截断误差的阶为 $O\left(h^{2}\right)$

将 $t\in [t_n,t_{n + 1}]$ 表成 $t = t_n + \tau h,0\leqslant \tau \leqslant 1.$ 由线性插值的余项公式，我们有

$$
\begin{array}{l} f (t, u (t)) = u ^ {\prime} (t) = u ^ {\prime} \left(t _ {n} + \tau h\right) \\ = u ^ {\prime} \left(t _ {n}\right) + \tau \left[ u ^ {\prime} \left(t _ {n + 1}\right) - u ^ {\prime} \left(t _ {n}\right) \right] + \frac {h ^ {2}}{2} \tau (\tau - 1) u ^ {\prime \prime \prime} \left(t _ {n} + \theta h\right) \quad (0 \leqslant \theta \leqslant 1). \\ \end{array}
$$

于是

$$
\begin{array}{l} \int_ {t _ {n}} ^ {t _ {n + 1}} f (t, u (t)) \mathrm {d} t = \int_ {0} ^ {1} \left[ u ^ {\prime} \left(t _ {n}\right) + \tau \left(u ^ {\prime} \left(t _ {n + 1}\right) - u ^ {\prime} \left(t _ {n}\right)\right) \right] h \mathrm {d} \tau + \\ \frac {h ^ {3}}{2} \int_ {0} ^ {1} \tau (\tau - 1) u ^ {\prime \prime \prime} (t _ {n} + \theta h) d \tau \\ = \frac {h}{2} \left[ u ^ {\prime} \left(t _ {n}\right) + u ^ {\prime} \left(t _ {n + 1}\right) \right] - \frac {h ^ {3}}{1 2} u ^ {\prime \prime \prime} (\zeta), \quad \zeta \in \left(t _ {n}, t _ {n + 1}\right). \\ \end{array}
$$

足见改进 Euler 法的局部截断误差为

$$
R _ {n} ^ {(1)} = - \frac {h ^ {3}}{1 2} u ^ {\prime \prime \prime} (\zeta), \tag {1.11}
$$

其阶为 $O(h^3)$ ，比Euler法高一阶

当然我们更关心的是近似解的误差，即

$$
e _ {n} = u \left(t _ {n}\right) - u _ {n},
$$

称为整体误差. 将 (1.4) 和 (1.8) 相减, 知 $e_n$ 满足误差方程:

$$
e _ {n + 1} = e _ {n} + h [ f (t _ {n}, u (t _ {n})) - f (t _ {n}, u _ {n}) ] + R _ {n}. \tag {1.12}
$$

因 $f(t,u)$ 关于 $\pmb{u}$ 满足Lipschitz条件(1.2)，故

$$
\begin{array}{l} \left| e _ {n + 1} \right| \leqslant \left| e _ {n} \right| + L h \left| e _ {n} \right| + R \\ = (1 + L h) \left| e _ {n} \right| + R, \\ \end{array}
$$

其中 $R = \max_{n}|R_{n}|$ 以此递推，得

$$
\begin{array}{l} \left| e _ {n} \right| \leqslant (1 + L h) \left| e _ {n - 1} \right| + R \leqslant (1 + L h) ^ {2} \left| e _ {n - 2} \right| + (1 + L h) R + R \\ \leqslant \dots \leqslant (1 + L h) ^ {n} | e _ {0} | + R \sum_ {j = 0} ^ {n - 1} (1 + L h) ^ {j} \\ = (1 + L h) ^ {n} \left| e _ {0} \right| + \frac {R}{L h} [ (1 + L h) ^ {n} - 1 ]. \\ \end{array}
$$

注意 $t_n = t_0 + nh \leqslant T, n = \frac{t_n - t_0}{h}$ , 于是

$$
\left| e _ {n} \right| \leqslant \mathrm {e} ^ {L (T - t _ {0})} \left| e _ {0} \right| + \frac {R}{L h} \left(\mathrm {e} ^ {L (T - t _ {0})} - 1\right), \quad n = 1, 2, \dots , N. \tag {1.13}
$$

右端依赖初始误差 $e_0$ 和局部截断误差的界 $R$ . 对 Euler 法, 可取 $R = Ch^2$ ( $C$ 是与 $n$ 无关的常数). 若 $e_0 = 0$ (取 $u_0 = u(t_0)$ ), 则

$$
\left| e _ {n} \right| \leqslant C L ^ {- 1} \mathrm {e} ^ {L (T - t _ {0})} h. \tag {1.14}
$$

所以 $e_n = O(h)$ ，比局部截断误差低一阶.用同样方法可以证明改进的Euler法的整体误差的阶为 $O\left(h^{2}\right)$ ，也比局部截断误差低一阶.

在实际计算中, 初值 $u_{0}$ 往往不能精确给出 (例如, 测量误差, 舍入误差等), 其误差将依次传递下去. 如果传递误差能够被控制, 精确说来, 传递误差连续依赖初始误差, 则称算法稳定; 否则就称不稳定. 显然不稳定的算法是不能用的. 我们考察 Euler 法. 设从初值 $u_{0}$ 和 $v_{0}$ 算出的节点值分别为 $\{u_{n}\}$ 和 $\{v_{n}\}$ , 则

$$
\begin{array}{l} u _ {n} = u _ {n - 1} + h f \left(t _ {n - 1}, u _ {n - 1}\right), \\ v _ {n} = v _ {n - 1} + h f \left(t _ {n - 1}, v _ {n - 1}\right), \quad n = 1, 2, \dots , N. \\ \end{array}
$$

两式相减并令 $e_n = u_n - v_n$ ，得

$$
e _ {n} = e _ {n - 1} + h \left[ f \left(t _ {n - 1}, u _ {n - 1}\right) - f \left(t _ {n - 1}, v _ {n - 1}\right) \right],
$$

从而

$$
\left| e _ {n} \right| \leqslant \left| e _ {n - 1} \right| + L h \left| e _ {n - 1} \right| = (1 + L h) \left| e _ {n - 1} \right|
$$

$$
\begin{array}{l} \leqslant \dots \leqslant (1 + L h) ^ {n} | e _ {0} | \\ \leqslant \mathrm {e} ^ {L T} \left| e _ {0} \right| \quad (\text {因} n h \leqslant T). \\ \end{array}
$$

这说明 $e_n$ 连续依赖初始误差 $e_0$ , 即 Euler 法稳定. 同样可证改进的 Euler 法也稳定.

# 1.1.3 线性差分方程

设 $a_0(n), a_1(n), \dots, a_k(n)$ 和 $b_n(n = 0,1,\dots)$ 满足 $a_0(n) \neq 0, a_k(n) \neq 0$ . 称序列 $\{u_n\}$ 满足的方程

$$
a _ {k} (n) u _ {n + k} + a _ {k - 1} (n) u _ {n + k - 1} + \dots + a _ {0} (n) u _ {n} = b _ {n}, \quad n = 0, 1, 2, \dots \tag {1.15}
$$

为 $k$ 阶线性差分方程，序列 $\{u_n\}$ 是差分方程的解。当右端 $b_{n} = 0(n = 0,1,\dots)$ 时，称为齐方程。为确定差分解，需给定 $k$ 个初值 $u_0,u_1,\dots ,u_{k - 1}$

记 $\Delta_{+}u_{n} = u_{n + 1} - u_{n}$ ，称为向前差分，则 $u_{n + 1} = u_n + \Delta_+u_n$ ，即 $u_{n + 1}$ 可用 $u_{n}$ 的一阶差分表示.又二阶差分

$$
\Delta_ {+} ^ {2} u _ {n} = \Delta_ {+} u _ {n + 1} - \Delta_ {+} u _ {n} = u _ {n + 2} - u _ {n + 1} - \Delta_ {+} u _ {n} = u _ {n + 2} - u _ {n} - 2 \Delta_ {+} u _ {n},
$$

故 $u_{n + 2} = u_n + 2\Delta_+u_n + \Delta_+^2 u_n$ ，即 $u_{n + 2}$ 可用 $u_{n}$ 的一阶和二阶差分表示.依次类推，可知 $u_{n + j}$ 能用 $u_{n}$ 的一阶、二阶直至 $j$ 阶差分表示.所以差分方程(1.15）的最高阶为 $k$ ： $k$ 阶线性差分方程是 $k$ 阶线性常微分方程的离散模拟，二者有许多平行的基本性质.例如：

(1) 齐方程的解具有可加性和齐次性. 若 $\{u_n\}$ 和 $\{v_n\}$ 都是齐方程的解, $\alpha$ 和 $\beta$ 是任意常数, 则 $\{\alpha u_n + \beta v_n\}$ 也是它的解.  
(2) $k$ 阶齐方程存在 $k$ 个线性无关的解. $k$ 个解 $\left\{u_n^{(j)}\right\} (j = 0,1,\dots ,k - 1)$ 称为线性无关，是指方程

$$
\sum_ {j = 0} ^ {k - 1} c _ {j} u _ {n} ^ {(j)} = 0, \quad n = 0, 1, \dots
$$

仅当 $c_{0} = c_{1} = \dots = c_{k - 1} = 0$ 时成立.由于任一解 $\left\{u_n^{(j)}\right\} (n = 0,1,\dots)$ 可表为初值 $u_0^{(j)},u_1^{(j)},\dots ,u_{k - 1}^{(j)}$ 的一次组合，所以 $k$ 个解 $\left\{u_n^{(j)}\right\} (n = 0,1,\dots ,k - 1)$ 线性无关的充要条件是初始向量 $\left(u_0^{(j)},u_1^{(j)},\dots ,u_{k - 1}^{(j)}\right)^{\mathrm{T}}$ （T表示转置， $j = 0,1,\dots ,k - 1)$ 线性无关，即行列式

$$
\left| \begin{array}{c c c c} u _ {0} ^ {(0)} & u _ {0} ^ {(1)} & \dots & u _ {0} ^ {(k - 1)} \\ u _ {1} ^ {(0)} & u _ {1} ^ {(1)} & \dots & u _ {1} ^ {(k - 1)} \\ \vdots & \vdots & & \vdots \\ u _ {k - 1} ^ {(0)} & u _ {k - 1} ^ {(1)} & \dots & u _ {k - 1} ^ {(k - 1)} \end{array} \right| \neq 0.
$$

由此进一步推出， $k$ 阶齐方程恰有 $k$ 个线性无关的解，且任一解可表示成这些解的线性组合.

(3) 非齐方程的通解等于齐方程的通解与非齐方程一特解之和.

今考虑常系数差分方程：

$$
\sum_ {j = 0} ^ {k} a _ {j} u _ {n + j} = b _ {n}, \quad n = 0, 1, \dots , \tag {1.16}
$$

其齐方程为

$$
\sum_ {j = 0} ^ {k} a _ {j} u _ {n + j} = 0, \quad n = 0, 1, \dots . \tag {1.17}
$$

考虑齐方程形如 $u_{n} = \zeta^{n}$ （ $\zeta$ 待定）的解，以之代到(1.17)，知 $\zeta$ 应满足

$$
a _ {k} \zeta^ {k} + a _ {k - 1} \zeta^ {k - 1} + \dots + a _ {1} \zeta + a _ {0} = 0,
$$

即 $\zeta$ 应是代数方程

$$
a _ {k} \lambda^ {k} + a _ {k - 1} \lambda^ {k - 1} + \dots + a _ {1} \lambda + a _ {0} = 0 \tag {1.18}
$$

的根. 反之, 若 $\zeta$ 是 (1.18) 的任一根, 则 $u_{n} = \zeta^{n}$ 必为 (1.17) 的解. 分几种情况:

(i) 方程 (1.18) 有 $k$ 个互异的实根 $\zeta_1, \zeta_2, \dots, \zeta_k$ , 则 $\zeta_1^n, \zeta_2^n, \dots, \zeta_k^n$ 是差分方程 (1.17) 的 $k$ 个线性无关解, 通解为

$$
u _ {n} = \sum_ {j = 1} ^ {k} c _ {j} \zeta_ {j} ^ {n}, \quad n = 0, 1, \dots .
$$

(ii) 方程 (1.18) 有 $m$ 个互异实根 $\zeta_1, \zeta_2, \dots, \zeta_m, \zeta_j$ 的重数是 $r_j (j = 1, 2, \dots, m)$ , $r_1 + r_2 + \dots + r_m = k$ , 则

$$
\zeta_ {j} ^ {n}, n \zeta_ {j} ^ {n}, \dots , n ^ {r _ {j} - 1} \zeta_ {j} ^ {n}, \quad n = 0, 1, \dots
$$

是 (1.17) 的 $r_j$ 个线性无关解. (1.17) 的通解形如

$$
u _ {n} = \sum_ {j = 1} ^ {m} \sum_ {l = 1} ^ {r _ {j}} c _ {j l} n ^ {l - 1} \zeta_ {j} ^ {n}, \quad n = 0, 1, \dots . \tag {1.19}
$$

(iii) 若 (1.18) 有复根 $\zeta_{j}$ , 则其共轭 $\overline{\zeta_{j}}$ 也是根. 令

$$
\zeta_ {j} = \rho \mathrm {e} ^ {\mathrm {i} \theta} = \rho (\cos \theta + \mathrm {i} \sin \theta) (i = \sqrt {- 1}),
$$

则

$$
\overline {{\zeta_ {j}}} = \overline {{\rho \mathrm {e} ^ {\mathrm {i} \theta}}} = \rho (\cos \theta - \mathrm {i} \sin \theta).
$$

此时可用两个线性无关的实解

$$
\rho^ {n} \cos n \theta , \rho^ {n} \sin n \theta
$$

替换 $\zeta_j^n$ 和 $\overline{\zeta_j^n}$

现在给出非齐方程(1.16)的通解表达式. 引进 $k$ 维向量 $\pmb{U}_n = (u_{n + k - 1}, u_{n + k - 2}, \dots, u_n)^{\mathrm{T}}$ 和 $k \times k$ 矩阵:

$$
C = \left[ \begin{array}{c c c c c} - a _ {k} ^ {- 1} a _ {k - 1} & - a _ {k} ^ {- 1} a _ {k - 2} & \dots & - a _ {k} ^ {- 1} a _ {1} & - a _ {k} ^ {- 1} a _ {0} \\ 1 & 0 & \dots & 0 & 0 \\ 0 & 1 & \dots & 0 & 0 \\ \vdots & \vdots & & \vdots & \vdots \\ 0 & 0 & \dots & 1 & 0 \end{array} \right]. \tag {1.20}
$$

将(1.16)改写成：

$$
u _ {n + k} = - a _ {k} ^ {- 1} \left(a _ {k - 1} u _ {n + k - 1} + a _ {k - 2} u _ {n + k - 2} + \dots + a _ {0} u _ {n}\right) + a _ {k} ^ {- 1} b _ {n},
$$

进一步写成向量形式：

$$
\boldsymbol {U} _ {n + 1} = \boldsymbol {C} \boldsymbol {U} _ {n} + \boldsymbol {b} _ {n}, \quad \boldsymbol {b} _ {n} = \left(a _ {k} ^ {- 1} b _ {n}, 0, \dots , 0\right) ^ {\mathrm {T}}. \tag {1.21}
$$

以此递推，即得通解

$$
\boldsymbol {U} _ {n} = \boldsymbol {C} ^ {n} \boldsymbol {U} _ {0} + \sum_ {l = 0} ^ {n - 1} \boldsymbol {C} ^ {l} \boldsymbol {b} _ {n - l - 1}, \quad n = 1, 2, \dots , \tag {1.22}
$$

其中第一项是齐方程的通解，第二项是非齐方程的特解（初值 $u_0 = u_1 = \dots = u_{k - 1} = 0$ ）

直接展开行列式 $|C - \lambda I_k|(I_k$ 是 $k$ 阶单位矩阵)，即知（1.18）的左端就是 $C$ 的特征多项式.设 $\lambda_{j}$ 是 $C$ 的特征值(方程(1.18)的根)， $x_{j} = (d_{k - 1},d_{k - 2},\dots ,d_{0})^{\mathrm{T}}$ 是相应的特征向量，则

$$
\begin{array}{l} - a _ {k} ^ {- 1} \left(a _ {k - 1} d _ {k - 1} + a _ {k - 2} d _ {k - 2} + \dots + a _ {0} d _ {0}\right) = \lambda_ {j} d _ {k - 1}, \\ d _ {k - 1} = \lambda_ {j} d _ {k - 2}, \\ \dots \\ d _ {1} = \lambda_ {j} d _ {0}. \\ \end{array}
$$

由此得

$$
\begin{array}{l} d _ {0} = d _ {0}, \quad d _ {1} = \lambda_ {j} d _ {0}, \quad d _ {2} = \lambda_ {j} ^ {2} d _ {0}, \quad \dots , \quad d _ {k - 1} = \lambda_ {j} ^ {k - 1} d _ {0}, \\ x _ {j} = d _ {0} \left(\lambda_ {j} ^ {k - 1}, \lambda_ {j} ^ {k - 2}, \dots , \lambda_ {j}, 1\right) ^ {\mathrm {T}}. \\ \end{array}
$$

可见任一特征值的特征空间的维数都是1，因此只有单特征值的初等因子的次数才是1.用相似变换 $S$ 将 $C$ 化成Jordan标准形：

$$
\boldsymbol {C} = \boldsymbol {S J S} ^ {- 1}.
$$

与单特征值相应的Jordan块为 $\lambda$ ，与重特征值相应的Jordan块为

$$
J _ {r} = \left[ \begin{array}{c c c c c} {\lambda} & {1} & & & {0} \\ & & & {\ddots} & \\ & & {\ddots} & & \\ {0} & & & {1} & \\ & & & {\lambda} \end{array} \right] (r: \lambda \text {的 重 数}).
$$

因 $C^n = SJ^n S^{-1},J^n$ 也是分块矩阵，每一分块形如 $(\lambda)^{n} = (\lambda^{n})$ （ $\lambda$ 是单特征值)或（是重特征值）

$$
J _ {r} ^ {n} = \left[ \begin{array}{c c c c} \lambda^ {n} & n \lambda^ {n - 1} & \dots & n \lambda^ {n - r - 1} \\ & \lambda^ {n} & \ddots & \vdots \\ & \mathbf {0} & \ddots & n \lambda^ {n - 1} \\ & & & \lambda^ {n} \end{array} \right] \quad (n \geqslant r).
$$

引理1.1 (i) 矩阵族 $\{C^n\} (n = 1,2,\dots)$ 有界的充要条件是：方程(1.18)的所有根在单位圆内 $(|\lambda |\leqslant 1)$ ，而位于单位圆周上的都是单根。

(ii) 矩阵族 $\{C^n\}$ 当 $n \to \infty$ 时有极限的充要条件是: 方程 (1.18) 的所有根在单位圆内, 而位于单位圆周上的根等于 1.

(iii) 矩阵族 $\{C^n\}$ 当 $n \to \infty$ 时趋于零矩阵的充要条件是: 方程 (1.18) 的所有根在单位圆内部 $(|\lambda| < 1)$ .

# 1.1.4 Gronwall 不等式

在做解的先验估计时经常要用Gronwall不等式（也称Bellman不等式). 先介绍连续形式的Gronwall不等式.

引理1.2 设连续函数 $\eta (t)(a\leqslant t\leqslant b)$ 满足

$$
| \eta (t) | \leqslant \beta + \alpha \int_ {a} ^ {t} | \eta (\tau) | \mathrm {d} \tau , \quad a \leqslant t \leqslant b, \tag {1.23}
$$

其中 $\alpha, \beta$ 为非负常数，则

$$
| \eta (t) | \leqslant \beta \mathrm {e} ^ {\alpha (t - a)}, \quad a \leqslant t \leqslant b. \tag {1.24}
$$

证明 先设 $\beta > 0$ . 令

$$
\zeta (t) = \beta + \alpha \int_ {0} ^ {t} | \eta (\tau) | \mathrm {d} \tau ,
$$

则由 (1.23) 得

$$
\frac {\mathrm {d} \zeta (t)}{\mathrm {d} t} \leqslant \alpha \zeta (t).
$$

显然 $\zeta (t) > 0$ ，故

$$
\frac {\zeta^ {\prime} (t)}{\zeta (t)} \leqslant \alpha .
$$

于 $[a,t]$ 上积分，得

$$
\ln \frac {\zeta (t)}{\beta} \leqslant \alpha (t - a).
$$

利用 (1.23) 即得 (1.24).

今设 $\beta = 0, \forall \delta > 0$ ，由（1.23）得

$$
| \eta (t) | \leqslant \delta + \alpha \int_ {a} ^ {t} | \eta (\tau) | d \tau .
$$

在上面得到的不等式中令 $\delta \to 0$ 便知 $\eta (t)\equiv 0,$ 故(1.24)仍成立

□

现在介绍离散形式的Gronwall不等式

引理1.3 设 $\alpha, \beta \geqslant 0$ 是任意常数，序列 $\{\eta_n\}$ 满足

$$
\left| \eta_ {n} \right| \leqslant \beta + \alpha h \sum_ {j = 0} ^ {n - 1} \eta_ {j}, \quad n = k, k + 1, \dots , n h \leqslant T, \tag {1.25}
$$

其中 $h > 0$ 是步长，则

$$
\left| \eta_ {n} \right| \leqslant \mathrm {e} ^ {\alpha T} (\beta + \alpha k h M _ {0}), \quad n \geqslant k, \quad n h \leqslant T, \tag {1.26}
$$

其中 $M_0 = \max \{| \eta_0 |, |\eta_1|, \dots, |\eta_{k-1}|\}$ .

证明 令 $\zeta_{n} = \beta +\alpha h\sum_{j = 0}^{n - 1}|\eta_{j}|$ ，则(1.25）相当于

$$
\zeta_ {n} - \zeta_ {n - 1} \leqslant \alpha h \zeta_ {n - 1},
$$

从而

$$
\zeta_ {n} \leqslant (1 + \alpha h) \zeta_ {n - 1} \leqslant (1 + \alpha h) ^ {2} \zeta_ {n - 2} \leqslant \dots \leqslant (1 + \alpha h) ^ {n - k} \zeta_ {k} \leqslant \mathrm {e} ^ {\alpha T} (\beta + \alpha k h M _ {0}),
$$

于是由 (1.25) 得

$$
\left| \eta_ {n} \right| \leqslant \zeta_ {n} \leqslant \mathrm {e} ^ {\alpha T} (\beta + \alpha k h M _ {0}).
$$

# 1.1.5 习题

1. 用 Euler 法和改进的 Euler 法求 $u' = -5u$ ( $0 \leqslant t \leqslant 1$ ), $u(0) = 1$ 的数值解, 步长 $h = 0.1, 0.05$ ; 并比较两个算法的精度.

2. 求差分方程 $u_{n+2} - 2\mu u_{n+1} + \mu u_n = 1 (n = 0,1,\dots)$ 的通解， $0 < \mu < 1$ 。证明 $u_n \to \frac{1}{1 - \mu}, n \to \infty$ 。  
3. 将 $u'' = -u$ ( $0 \leqslant t \leqslant 1$ ), $u(0) = 0$ , $u'(0) = 1$ 化为一阶方程组, 并用 Euler 法和改进的 Euler 法求解, 步长 $h = 0.1, 0.05$ ; 并比较两个算法的精度.

# 1.2 线性多步法

用Euler法计算节点 $t_n = t_0 + nh(t_0 = 0)$ 的近似值 $u_{n}$ 只用到前一节点的值 $u_{n - 1}$ 所以从初值 $u_{0}$ 出发可算出以后各节点的值，这样的方法称为单步法.为了提高解的精度，需要构造线性多步法，其一般形式为

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u _ {n + j} = h \sum_ {j = 0} ^ {k} \beta_ {j} f _ {n + j}, \tag {1.27}
$$

其中 $f_{n + j} = f(t_{n + j},u_{n + j}),\alpha_{j}$ 和 $\beta_{j}$ 是常数，且 $\alpha_{k}\neq 0,\alpha_{0}$ 和 $\beta_0$ 不同时为0.按(1.27)计算 $u_{n + k}$ 时要用到前面 $k$ 个节点的值 $u_{n},u_{n + 1},\dots ,u_{n + k - 1}$ ，因此称（1.27）为多步法或 $k$ 步法.又因为方程(1.27)关于 $f_{n + j}$ 是线性的，所以称为线性多步法.为使多步法的计算能够进行，除给定的初值 $u_{0}$ 外，还要知道附加初值 $u_{1},u_{2},\dots ,u_{k - 1}$ ，这可用其他方法计算，后面还要介绍.由于多步法每计算一步用到的信息更多，因此可望造出精度更高的算法.若 $\beta_{k} = 0$ ，则称方法(1.27)是显式的；若 $\beta_{k}\neq 0$ ，则称方法(1.27）是隐式的.

# 1.2.1 数值积分法

将方程 $u^{\prime} = f(t,u)$ 写成积分形式，比如在 $[t_n,t_{n + 1}]$ 上积分，得

$$
u \left(t _ {n + 1}\right) = u \left(t _ {n}\right) + \int_ {t _ {n}} ^ {t _ {n + 1}} f (t, u (t)) \mathrm {d} t. \tag {1.28}
$$

适当取 $k + 1$ 个节点，以被积函数 $f(t,u(t))$ 的 $k$ 次Lagrange插值多项式 $L_{n,k}(t)$ 近似代替 $f(t,u(t))$ ，就可得到形如（1.27）的线性多步法．插值节点的不同取法就会得到不同的多步法.

(1) Adams外插法也称Adams-Bashforth法, 这是一种显式多步法. 取 $t_n, t_{n-1}, \cdots, t_{n-k}$ 为节点, 构造 $f$ 的Langrange插值多项式 $L_{n,k}(t)$ , 则

$$
f (t, u (t)) = L _ {n, k} (t) + r _ {n, k} (t), \tag {1.29}
$$

其中 $r_{n,k}(t)$ 是插值余项. 代到 (1.28), 得

$$
u \left(t _ {n + 1}\right) = u \left(t _ {n}\right) + \int_ {t _ {n}} ^ {t _ {n + 1}} L _ {n, k} (t) \mathrm {d} t + \int_ {t _ {n}} ^ {t _ {n + 1}} r _ {n, k} (t) \mathrm {d} t. \tag {1.30}
$$

舍去余项

$$
R _ {n, k} = \int_ {t _ {n}} ^ {t _ {n + 1}} r _ {n, k} (t) \mathrm {d} t, \tag {1.31}
$$

并用 $u_{j}$ 代替 $u(t_j)$ ，即得

$$
u _ {n + 1} = u _ {n} + \int_ {t _ {n}} ^ {t _ {n + 1}} L _ {n, k} (t) \mathrm {d} t. \tag {1.32}
$$

像 Euler 法一样, 称 $R_{n,k}$ 为局部截断误差.

现在给出 (1.32) 的具体形式. 因为插值节点等距, 被插值点 $t \in [t_n, t_{n+1}]$ 靠近最后一个节点 $t_n$ , 所以将 $L_{n,k}(t)$ 表示成牛顿向后插值公式更方便. 记

$$
t = t _ {n} + \tau h, \quad \tau \in [ 0, 1 ],
$$

则牛顿向后插值公式为

$$
\begin{array}{l} L _ {n, k} (t) = L _ {n, k} \left(t _ {n} + \tau h\right) \\ = f _ {n} + \frac {\tau}{1 !} \Delta_ {+} f _ {n - 1} + \frac {\tau (\tau + 1)}{2 !} \Delta_ {+} ^ {2} f _ {n - 2} + \dots + \frac {\tau (\tau + 1) \cdots (\tau + k - 1)}{k !} \Delta_ {+} ^ {k} f _ {n + k}, \tag {1.33} \\ \end{array}
$$

式中 $\Delta_{+}^{j}$ 表示 $j$ 阶向前差分， $f_{n-j} = f(t_{n-j}, u_{n-j})$ 。引进记号

$$
\binom {s} {j} = \frac {s (s - 1) \cdots (s - j + 1)}{j !}, \quad \binom {s} {0} = 1, \tag {1.34}
$$

则

$$
L _ {n, k} (t) = \sum_ {j = 0} ^ {k} (- 1) ^ {j} \binom {- \tau} {j} \Delta_ {+} ^ {j} f _ {n - j}. \tag {1.35}
$$

以之代到 (1.32), 得

$$
u _ {n + 1} = u _ {n} + h \sum_ {j = 0} ^ {k} a _ {j} \Delta_ {+} ^ {j} f _ {n - j}, \tag {1.36}
$$

其中

$$
a _ {j} = (- 1) ^ {j} \int_ {0} ^ {1} \binom {- \tau} {j} d \tau , \quad j = 0, 1, \dots , k. \tag {1.37}
$$

令 $\eta_{j} = hf_{j}$ ，则(1.36)还可写成

$$
u _ {n + 1} = u _ {n} + \sum_ {j = 0} ^ {k} a _ {j} \Delta_ {+} ^ {j} \eta_ {n - j}. \tag {1.38}
$$

这就是 Adams 外插公式. 显然当 $k = 0$ 时就是 Euler 法

为了计算 $a_{j}$ ，我们给出联系这些系数的递推公式.定义 $\{a_j\}_{0}^{\infty}$ 的生成函数

$$
G (t) = \sum_ {j = 0} ^ {\infty} a _ {j} t ^ {j},
$$

其中 $a_{j}$ 同(1.37)， $j = 1,2,\dots$ 因为 $a_{j}$ 的界不超过1，故级数当 $|t| < 1$ 时收敛.将(1.37)代到上式，则

$$
\begin{array}{l} G (t) = \sum_ {j = 0} ^ {\infty} (- 1) ^ {j} \int_ {0} ^ {1} \binom {- \tau} {j} t ^ {j} d \tau = \int_ {0} ^ {1} \sum_ {j = 0} ^ {\infty} (- t) ^ {j} \binom {- \tau} {j} d \tau \\ = \int_ {0} ^ {1} (1 - t) ^ {- \tau} d \tau = - \frac {t}{(1 - t) \ln (1 - t)}, \\ \end{array}
$$

从而

$$
- \frac {\ln (1 - t)}{t} G (t) = \frac {1}{1 - t}.
$$

两端展成幂级数，上式就是

$$
\left(1 + \frac {1}{2} t + \frac {1}{3} t ^ {2} + \dots\right) \left(a _ {0} + a _ {1} t + a _ {2} t ^ {2} + \dots\right) = 1 + t + t ^ {2} + \dots .
$$

比较 $t^n$ 的系数，则得递推式：

$$
a _ {n} + \frac {1}{2} a _ {n - 1} + \frac {1}{3} a _ {n - 2} + \dots + \frac {1}{n + 1} a _ {0} = 1, \quad n = 0, 1, \dots .
$$

由此可依次算出系数 $a_{j}$ 。表1.1给出前几个系数值

表 1.1 Adams 外插公式系数的值  

<table><tr><td>j</td><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td></tr><tr><td>aj</td><td>1</td><td>1/2</td><td>5/12</td><td>3/8</td><td>251/720</td><td>95/288</td><td>10987/60480</td></tr></table>

回想插值公式的余项为

$$
r _ {n, k} (t) = r _ {n, k} \left(t _ {n} + \tau h\right) = (- 1) ^ {k + 1} \binom {- \tau} {k + 1} h ^ {k + 1} u ^ {(k + 2)} (\bar {\xi}), \tag {1.39}
$$

其中 $t_{n - k}\leqslant \bar{\xi}\leqslant t_{n + 1}$ ，则知

$$
\begin{array}{l} R _ {n, k} = h ^ {k + 2} \int_ {0} ^ {1} (- 1) ^ {k + 1} \left( \begin{array}{c} - \tau \\ k + 1 \end{array} \right) u ^ {(k + 2)} (\bar {\xi}) \mathrm {d} \tau \\ = a _ {k + 1} h ^ {k + 2} u ^ {(k + 2)} (\xi), t _ {n - k} \leqslant \xi \leqslant t _ {n + 1}. \tag {1.40} \\ \end{array}
$$

这里用到了积分第二中值公式. 由 (1.40) 知道 Adams 外插法 (1.32) 或 (1.36) 的局部截

断误差的阶为 $O\left(h^{k + 2}\right)$

实际计算时，常常将(1.36)右端的差分表为 $f_{n - j}$ 的线性组合.为此利用差分公式

$$
\Delta_ {+} ^ {j} f _ {n - j} = \sum_ {l = 0} ^ {j} (- 1) ^ {l} \binom {j} {l} f _ {n - l}, \tag {1.41}
$$

则(1.36)化为

$$
u _ {n + 1} = u _ {n} + h \sum_ {l = 0} ^ {k} b _ {k l} f _ {n - l}, \tag {1.42}
$$

其中

$$
b _ {k l} = (- 1) ^ {l} \sum_ {j = l} ^ {k} a _ {j} \binom {j} {l}. \tag {1.43}
$$

利用表1.1，可给出 $b_{kl}$ 如表1.2

表 1.2 系数 ${b}_{kl}$ 的值  

<table><tr><td>l</td><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td></tr><tr><td>b0l</td><td>1</td><td></td><td></td><td></td><td></td><td></td></tr><tr><td>2b1l</td><td>3</td><td>-1</td><td></td><td></td><td></td><td></td></tr><tr><td>12b2l</td><td>23</td><td>-16</td><td>5</td><td></td><td></td><td></td></tr><tr><td>24b3l</td><td>55</td><td>-59</td><td>37</td><td>-9</td><td></td><td></td></tr><tr><td>720b4l</td><td>1901</td><td>-2774</td><td>2616</td><td>-1274</td><td>251</td><td></td></tr><tr><td>1440b5l</td><td>4277</td><td>-7923</td><td>9982</td><td>-7298</td><td>2877</td><td>-475</td></tr></table>

例如 $k = 0,1,2,3$ 的外插公式分别为

$$
\begin{array}{l} k = 0: \quad u _ {n + 1} = u _ {n} + h f \left(t _ {n}, u _ {n}\right), \\ k = 1: \quad u _ {n + 1} = u _ {n} + \frac {h}{2} \left(3 f _ {n} - f _ {n - 1}\right), \\ k = 2: \quad u _ {n + 1} = u _ {n} + \frac {h}{1 2} \left(2 3 f _ {n} - 1 6 f _ {n - 1} + 5 f _ {n - 2}\right), \\ k = 3: \quad u _ {n + 1} = u _ {n} + \frac {h}{2 4} (5 5 f _ {n} - 5 9 f _ {n - 1} + 3 7 f _ {n - 2} - 9 f _ {n - 3}). \\ \end{array}
$$

(2) Adams 内插法 也称 Adams-Moulton 法, 这是一种隐式多步法. 现在取插值节点为 $t_{n-k}, t_{n-k+1}, \cdots, t_n, t_{n+1}$ (比外插法多取一点 $t_{n+1}$ ), 构造 $u'(t)$ 或 $f(t, u(t))$ 的 $k+1$ 次 Lagrange 插值多项式 $L_{n,k}^{(1)}(t)$ , 插值余项为 $r_{n,k}^{(1)}(t)$ , 则

$$
f = L _ {n, k} ^ {(1)} (t) + r _ {n, k} ^ {(1)} (t).
$$

以之代到 (1.28) 右端, 得

$$
u \left(t _ {n + 1}\right) = u \left(t _ {n}\right) + \int_ {t _ {n}} ^ {t _ {n + 1}} L _ {n, k} ^ {(1)} (t) \mathrm {d} t + \int_ {t _ {n}} ^ {t _ {n + 1}} r _ {n, k} ^ {(1)} (t) \mathrm {d} t. \tag {1.44}
$$

舍去余项

$$
R _ {n, k} ^ {(1)} = \int_ {t _ {n}} ^ {t _ {n + 1}} r _ {n, k} ^ {(1)} (t) \mathrm {d} t, \tag {1.45}
$$

并用 $u_{j}$ 代替 $u(t_{j})$ ，则得Adams内插法：

$$
u _ {n + 1} = u _ {n} + \int_ {t _ {n}} ^ {t _ {n + 1}} L _ {n, k} ^ {(1)} (t) \mathrm {d} t. \tag {1.46}
$$

当 $k = 0$ 时, 就是改进的 Euler 法. 余项 $R_{n,k}^{(1)}$ 是内插法的局部截断误差

现在将（1.46）具体化.仍用牛顿向后插值公式

$$
L _ {n, k} ^ {(1)} (t) = L _ {n, k} ^ {(t)} (t) = \sum_ {j = 0} ^ {k + 1} (- 1) ^ {j} \binom {- \tau} {j} \Delta_ {+} ^ {j} f _ {n - j + 1}, \tag {1.47}
$$

其中 $\tau \in [-1,0]$ ，而二项系数

$$
\begin{array}{l} \binom {- \tau} {j} = \frac {- \tau (- \tau - 1) \cdots (- \tau - j + 1)}{j !} \\ = (- 1) ^ {j} \frac {\tau (\tau + 1) \cdots (\tau + j - 1)}{j !}. \tag {1.48} \\ \end{array}
$$

将(1.47)代到(1.46)右端，则得Adams内插公式：

$$
u _ {n + 1} = u _ {n} + h \sum_ {j = 0} ^ {k + 1} a _ {j} ^ {*} \Delta_ {+} ^ {j} f _ {n - j + 1} = u _ {n} + \sum_ {j = 0} ^ {k + 1} a _ {j} ^ {*} \Delta_ {+} ^ {j} \eta_ {n - j + 1}, \tag {1.49}
$$

其中 $\eta_{n - j + 1} = hf_{n - j + 1} = hf(t_{n - j + 1},u_{n - j + 1})$ ，而

$$
a _ {j} ^ {*} = (- 1) ^ {j} \int_ {- 1} ^ {0} \binom {- \tau} {j} d \tau , \quad j = 0, 1, \dots , k + 1. \tag {1.50}
$$

用生成函数法可导出系数 $a_j^*$ 的递推公式. 表 1.3 给出 $a_j^*$ 的前几个值

表 1.3 系数 ${a}_{j}^{ * }$ 的值  

<table><tr><td>j</td><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td></tr><tr><td>aj*</td><td>1</td><td>-1/2</td><td>-1/12</td><td>-1/24</td><td>-19/720</td><td>-3/160</td><td>-863/60480</td></tr></table>

利用插值余项公式

$$
\begin{array}{l} r _ {n, k} ^ {(1)} (t) = r _ {n, k} ^ {(1)} \left(t _ {n + 1} + \tau h\right) \\ = (- 1) ^ {k} \binom {- \tau} {k + 2} h ^ {k + 2} u ^ {(k + 3)} (\bar {\xi}), \quad t _ {n - k} \leqslant \bar {\xi} \leqslant t _ {n + 1}, \tag {1.51} \\ \end{array}
$$

则得

$$
R _ {n, k} ^ {(1)} = a _ {k + 2} ^ {*} h ^ {k + 3} u ^ {(k + 3)} (\xi), \quad t _ {n - k} \leqslant \xi \leqslant t _ {n + 1}. \tag {1.52}
$$

这里用到了积分第二中值公式.由此可见，Adams内插法的局部截断误差的阶为 $O\left(h^{k + 3}\right)$ 利用差商公式(1.41)，可将(1.49）写成便于计算的形式：

$$
u _ {n + 1} = u _ {n} + h \sum_ {l = 0} ^ {k + 1} b _ {k + 1, l} ^ {*} f _ {n - l + 1}, \tag {1.53}
$$

其中

$$
b _ {k + 1, l} ^ {*} = (- 1) ^ {l} \sum_ {j = l} ^ {k + 1} a _ {j} ^ {*} \binom {j} {l}. \tag {1.54}
$$

表 1.4 列出 ${b}_{k + 1,l}^{ * }$ 的值.

表 1.4 系数 ${b}_{k + 1,l}^{ * }$ 的值  

<table><tr><td>l</td><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td></tr><tr><td>b0l*</td><td>1</td><td></td><td></td><td></td><td></td><td></td></tr><tr><td>2b1l*</td><td>1</td><td>1</td><td></td><td></td><td></td><td></td></tr><tr><td>12b2l*</td><td>5</td><td>8</td><td>-1</td><td></td><td></td><td></td></tr><tr><td>24b3l*</td><td>9</td><td>19</td><td>-5</td><td>1</td><td></td><td></td></tr><tr><td>720b4l*</td><td>251</td><td>646</td><td>-264</td><td>106</td><td>-19</td><td></td></tr><tr><td>1440b5l*</td><td>475</td><td>1427</td><td>-798</td><td>482</td><td>-173</td><td>27</td></tr></table>

例如 $k = 0,1,2,3$ 的内插公式分别为

$$
\begin{array}{l} k = 0: \quad u _ {n + 1} = u _ {n} + h f _ {n + 1}, \\ k = 1: \quad u _ {n + 1} = u _ {n} + \frac {h}{2} \left(f _ {n + 1} + f _ {n}\right), \\ k = 2: \quad u _ {n + 1} = u _ {n} + \frac {h}{1 2} \left(5 f _ {n + 1} + 8 f _ {n} - f _ {n - 1}\right), \\ k = 3: \quad u _ {n + 1} = u _ {n} + \frac {h}{2 4} (9 f _ {n + 1} + 1 9 f _ {n} - 5 f _ {n - 1} + f _ {n - 2}). \\ \end{array}
$$

Adams外插法和内插法有以下几点区别

第一, 从表1.1和表1.3(表1.2和表1.4)可知, 按绝对值系数 $a_{j}^{*}$ 比 $a_{j}$ 小 $(b_{k + 1,l}^{*}$ 比 $b_{kl}$ 小), 因此计算中内插法的舍入误差影响比外插法小.

第二，用外插法和内插法计算 $t_{n + 1}$ 处的值 $u_{n + 1}$ ，用到的已知量相同 $(k + 1$ 个值 $u_{n}, u_{n - 1}, \dots, u_{n - k})$ ，但内插法局部截断误差的阶为 $O(h^{k + 3})$ ，外插法局部截断误差的阶为 $O(h^{k + 2})$ ，前者比后者高一阶。所以为达到相同的误差阶，内插法比外插法可少用一个初始已知量。

第三, 外插法是显式, 计算 $u_{n+1}$ 是直接的. 内插法是隐式, 计算 $u_{n+1}$ 需要解方程

(1.53), 通常用如下迭代求解:

$$
u _ {n + 1} ^ {[ m + 1 ]} = u _ {n} + h b _ {k + 1, 0} ^ {*} f \left(t _ {n + 1}, u _ {n + 1} ^ {[ m ]}\right) + h \sum_ {l = 1} ^ {k + 1} b _ {k + 1, l} ^ {*} f _ {n - l + 1}, \quad m = 0, 1 \dots . \tag {1.55}
$$

当 $h$ 充分小时，可使

$$
\left| h b _ {k + 1, 0} ^ {*} \frac {\partial f (t _ {n + 1} , u _ {n + 1})}{\partial u _ {n + 1}} \right| \leqslant h b _ {k + 1, 0} ^ {*} L <   1,
$$

此时迭代（1.55）收敛. 初始近似可用外插法给出，即

$$
u _ {n + 1} ^ {[ 0 ]} = u _ {n} + h \sum_ {l = 0} ^ {k} b _ {k l} f _ {n - l}. \tag {1.56}
$$

由于这是好的近似, 所以收敛是很快的, 通常迭代2至3次就可收敛.

Adams外插法和内插法是这样得到的，先将初值问题改写成积分形式(1.28)，再用适当的数值积分离散(1.28).其实，也可将初值问题写成其他积分形式，例如：

$$
u (t _ {n + 2}) - u (t _ {n}) = \int_ {t _ {n}} ^ {t _ {n + 2}} f (t, u (t)) \mathrm {d} t,
$$

再用适当的数值积分代替右端积分，例如用Simpson公式，得到

$$
u _ {n + 2} - u _ {n} = \frac {h}{3} \left(f _ {n} + 4 f _ {n + 1} + f _ {n + 2}\right), \tag {1.57}
$$

这是线性二步法.

还应指出，用数值积分法只能构造一类特殊的多步法，其系数

$$
\alpha_ {k} = 1, \quad \alpha_ {k - m} = - 1, \quad \alpha_ {l} = 0 \quad (l \neq k - m, k).
$$

下面介绍更一般的待定系数法

# 1.2.2 待定系数法

令

$$
L [ u (t); h ] = \sum_ {j = 0} ^ {k} \left[ \alpha_ {j} u (t + j h) - h \beta_ {j} u ^ {\prime} (t + j h) \right]. \tag {1.58}
$$

设 $u(t)$ 是初值问题的解，将 $u(t + jh)$ 和 $u^{\prime}(t + jh)$ 在点 $t$ 用Taylor公式展开，按 $h$ 的同次幂合并同类项，得

$$
L [ u (t); h ] = c _ {0} u (t) + c _ {1} h u ^ {(1)} (t) + \dots + c _ {q} h ^ {q} u ^ {(q)} (t) + \dots , \tag {1.59}
$$

其中

$$
\left\{ \begin{array}{l} c _ {0} = \alpha_ {0} + \alpha_ {1} + \dots + \alpha_ {k}, \\ c _ {1} = \alpha_ {1} + 2 \alpha_ {2} + \dots + k \alpha_ {k} - (\beta_ {0} + \beta_ {1} + \dots + \beta_ {k}), \\ \dots \dots \dots \dots . \\ c _ {q} = \frac {1}{q !} (\alpha_ {1} + 2 ^ {q} \alpha_ {2} + \dots + k ^ {q} \alpha_ {k}) - \\ \frac {1}{(q - 1) !} \cdot (\beta_ {1} + 2 ^ {q - 1} \beta_ {2} + \dots + k ^ {q - 1} \beta_ {k}), \quad q = 2, 3, \dots . \end{array} \right. \tag {1.60}
$$

若 $u(t)$ 有 $p + 2$ 次连续微商，则可选取 $k$ （足够大）和 $\alpha_{j},\beta_{j}$ 使 $c_{0} = c_{1} = \dots = c_{p} = 0,$ 而 $c_{p + 1}\neq 0$ ，即选择 $\alpha_{j},\beta_{j}$ 满足

$$
\left\{ \begin{array}{l} \alpha_ {0} + \alpha_ {1} + \dots + \alpha_ {k} = 0, \\ \alpha_ {1} + 2 \alpha_ {2} + \dots + k \alpha_ {k} - (\beta_ {0} + \beta_ {1} + \dots + \beta_ {k}) = 0, \\ \dots \dots \dots \dots . \\ \frac {1}{p !} (\alpha_ {1} + 2 ^ {p} \alpha_ {2} + \dots + k ^ {p} \alpha_ {k}) - \frac {1}{(p - 1) !} \cdot (\beta_ {1} + 2 ^ {p - 1} \beta_ {2} + \dots + k ^ {p - 1} \beta_ {k}) = 0. \end{array} \right. \tag {1.61}
$$

此时

$$
L [ u (t); h ] = c _ {p + 1} h ^ {p + 1} u ^ {(p + 1)} (t) + O \left(h ^ {p + 2}\right). \tag {1.62}
$$

由于 $u^{\prime}(t) = f(t,u(t))$ ，则

$$
\sum_ {j = 0} ^ {k} \left[ \alpha_ {j} u (t _ {n} + j h) - h \beta_ {j} f (t _ {n} + j h, u (t _ {n} + j h)) \right] = R _ {n, k}, \tag {1.63a}
$$

$$
R _ {n, k} = c _ {p + 1} h ^ {p + 1} u ^ {(p + 1)} \left(t _ {n}\right) + O \left(h ^ {p + 2}\right). \tag {1.63b}
$$

略去余项 $R_{n,k}$ , 并用 $u_{n+j}$ 代替 $u(t_n + jh)$ , 用 $f_{n+j}$ 记 $f(t_{n+j}, u_{n+j})$ , 就得到线性多步法 (1.27), 其局部截断误差 $R_{n,k} = O(h^{p+1})$ . 往后将证明方法的整体误差的阶是 $O(h^p)$ , 所以称此法为 $p$ 阶 $k$ 步法. 显然 $p$ 的大小和 $k$ 有关.

因为多步法 (1.27) 可以差一个非零乘数, 所以不妨设 $\alpha_{k} = 1$ . 当 $\beta_{k} = 0$ 时, $u_{n+k}$ 可用 $u_{n+k-1}, u_{n+k-2}, \cdots, u_{n}$ 直接表示, 此为显方法. 反之, 当 $\beta_{k} \neq 0$ 时, 求 $u_{n+k}$ 需解一个方程 (一般用迭代法), 此为隐方法. 用待定系数法构造多步法的一个基本要求, 是选取 $\alpha_{j}, \beta_{j}$ 使局部截断误差的阶尽可能高.

作为待定系数法的一个应用，我们讨论一般的二步法.此时 $k = 2,\alpha_{2} = 1.$ 记 $\alpha_0 =$ $\alpha$ ，其余四个系数 $\alpha_{1},\beta_{0},\beta_{1},\beta_{2}$ 由 $c_{0} = c_{1} = c_{2} = c_{3} = 0$ 确定，即满足方程：

$$
\left\{ \begin{array}{l} c _ {0} = \alpha + \alpha_ {1} + 1 = 0, \\ c _ {1} = \alpha_ {1} + 2 - (\beta_ {0} + \beta_ {1} + \beta_ {2}) = 0, \\ c _ {2} = \frac {1}{2} (\alpha_ {1} + 4) - (\beta_ {1} + 2 \beta_ {2}) = 0, \\ c _ {3} = \frac {1}{6} (\alpha_ {1} + 8) - \frac {1}{2} (\beta_ {1} + 4 \beta_ {2}) = 0. \end{array} \right.
$$

解之得

$$
\alpha_ {1} = - (1 + \alpha), \quad \beta_ {0} = - \frac {1}{1 2} (1 + 5 \alpha),
$$

$$
\beta_ {1} = \frac {2}{3} (1 - \alpha), \quad \beta_ {2} = \frac {1}{1 2} (5 + \alpha).
$$

所以一般二步法为

$$
u _ {n + 2} - (1 + \alpha) u _ {n + 1} + \alpha u _ {n} = \frac {h}{1 2} [ (5 + \alpha) f _ {n + 2} + 8 (1 - \alpha) f _ {n + 1} - (1 + 5 \alpha) f _ {n} ]. \tag {1.64}
$$

由 (1.60) 还知道

$$
c _ {4} = \frac {1}{2 4} (\alpha_ {1} + 1 6) - \frac {1}{6} (\beta_ {1} + 8 \beta_ {2}) = - \frac {1}{2 4} (1 + \alpha),
$$

$$
c _ {5} = \frac {1}{1 2 0} (\alpha_ {1} + 3 2) - \frac {1}{2 4} (\beta_ {1} + 1 6 \beta_ {2}) = - \frac {1}{3 6 0} (1 7 + 1 3 \alpha).
$$

所以当 $\alpha \neq -1$ 时, $c_{4} \neq 0$ , 方法 (1.64) 是三阶二步法. 当 $\alpha = -1$ 时, $c_{4} = 0$ , 但 $c_{5} \neq 0$ , 方法 (1.64) 化为

$$
u _ {n + 2} = u _ {n} + \frac {h}{3} \left(f _ {n + 2} + 4 f _ {n + 1} + f _ {n}\right), \tag {1.65}
$$

这是四阶二步法, 是具有最高阶的二步法, 称为Milne法. 前面我们曾用Simpson公式导出这一算法(见(1.57)). 此外若取 $\alpha = 0$ , 则(1.64)为二步Adams内插法; 若取 $\alpha = -5$ , 则(1.64)是显方法.

# 1.2.3 预估-校正算法

将隐式 $k$ 步法 (1.27) 写成:

$$
u _ {n + k} + \sum_ {j = 0} ^ {k - 1} \alpha_ {j} u _ {n + j} = h \beta_ {k} f (t _ {n + k}, u _ {n + k}) + h \sum_ {j = 0} ^ {k - 1} \beta_ {j} f _ {n + j}, \tag {1.66}
$$

其中 $f_{n + j} = f(t_{n + j},u_{n + j})$ .若已求出 $u_{n + j},j = 0,1,\dots ,k - 1,$ 则(1.66)关于 $u_{n + k}$ 为非线性方程，通常用如下迭代法求解：

$$
u _ {n + k} ^ {[ m + 1 ]} + \sum_ {j = 0} ^ {k - 1} \alpha_ {j} u _ {n + j} = h \beta_ {k} f (t _ {n + k}, u _ {n + k} ^ {[ m ]}) + h \sum_ {j = 0} ^ {k - 1} \beta_ {j} f _ {n + j}, m = 0, 1, \dots , \tag {1.67}
$$

其中 $u_{n + k}^{[0]}$ 为给定的迭代初值. 显然若

$$
h <   \frac {1}{L} | \beta_ {k} |,
$$

$L$ 为 $f$ 关于 $u$ 的Lipschitz常数，初值 $u_{n + k}^{[0]}$ 又选择适当，则迭代(1.67)收敛

因隐式方法 (1.66) 每步的计算量取决于迭代 (1.67) 的次数, 所以选好初值 $u_{n+k}^{[0]}$ 非

常重要.最自然的一种方法是用显式多步法计算 $u_{n + k}^{[0]}$ ，比如

$$
u _ {n + k} ^ {[ 0 ]} + \sum_ {j = 0} ^ {k - 1} \alpha_ {j} ^ {*} u _ {n + j} = h \sum_ {j = 0} ^ {k - 1} \beta_ {j} ^ {*} f _ {n + j}. \tag {1.68}
$$

称 (1.68) 为预估算式 ( $P$ 算式), (1.67) 为校正算式 ( $C$ 算式), 统称 (1.67)—(1.68) 为预估-校正算法, 简称预-校算法或 $PC$ 算法 (Predictor-Corrector methods).

一个极端情形是允许 (1.67) 不断进行, 直至不等式 $\left|u_{n+k}^{[m+1]} - u_{n+k}^{[m]}\right| < \varepsilon$ 成立, 其中 $\varepsilon$ 是预先指定的容许误差. 由于这种算法对迭代次数不加限制, 花费在计算函数 $f$ 的工作量可能很大, 所以通常采用另一种限制迭代次数的算法. 假定校正次数 $M$ (即迭代次数) 固定, $P$ 表示预估算子, $C$ 是一次校正算子 (即迭代一次), $E$ 是计算 $f$ 一次的运算, 则预估一次校正 $M$ 次的算法可记为 $P(EC)^{M} = P(EC)(EC) \cdots (EC)$ , 计算过程如下:

$$
\begin{array}{l} P: \quad u _ {n + k} ^ {[ 0 ]} + \sum_ {j = 0} ^ {k - 1} \alpha_ {j} ^ {*} u _ {n + j} ^ {[ M ]} = h \sum_ {j = 0} ^ {k - 1} \beta_ {j} ^ {*} f _ {n + j} ^ {[ M - 1 ]}, \\ E: \quad f _ {n + k} ^ {[ m ]} = f \left(t _ {n + k}, u _ {n + k} ^ {[ m ]}\right), \tag {1.69} \\ C: \quad u _ {n + k} ^ {[ m + 1 ]} + \sum_ {j = 0} ^ {k - 1} \alpha_ {j} u _ {n + j} ^ {[ M ]} = h \beta_ {k} f _ {n + k} ^ {[ m ]} + h \sum_ {j = 0} ^ {k - 1} \beta_ {j} f _ {n + j} ^ {[ M - 1 ]}, \\ \end{array}
$$

其中 $m = 0,1,\dots ,M - 1.$ 按这一预-校格式计算结束时，得到的数据是 $u_{n + k}^{[M]}$ 和 $f_{n + k}^{[M - 1]} =$ $f\left(t_{n + k},u_{n + k}^{[M - 1]}\right)$ ，为下一步 $(t = t_{n + 1 + k})$ 预估计算所用.显然 $u_{n + k}^{[M]}$ 比 $u_{n + k}^{[M - 1]}$ 更接近 $u_{n + k}$ ，因此还可以设计一种算法，每一步算出 $u_{n + k}^{[M]}$ 后，利用它将 $f_{n + k}^{[M]} = f\left(t_{n + k},u_{n + k}^{[M]}\right)$ 算出，供下一步预估时使用.这种预-校算法记为 $P(EC)^{M}E$ ，计算过程如下：

$$
\begin{array}{l} P: u _ {n + k} ^ {[ 0 ]} + \sum_ {j = 0} ^ {k - 1} \alpha_ {j} ^ {*} u _ {n + j} ^ {[ M ]} = h \sum_ {j = 0} ^ {k - 1} \beta_ {j} ^ {*} f _ {n + j} ^ {[ M ]}, \\ E: f _ {n + k} ^ {[ m ]} = f \left(t _ {n + k}, u _ {n + k} ^ {[ m ]}\right), \\ C: u _ {n + k} ^ {[ m + 1 ]} + \sum_ {j = 0} ^ {k - 1} \alpha_ {j} u _ {n + j} ^ {[ M ]} = h \beta_ {k} f _ {n + k} ^ {[ m ]} + h \sum_ {j = 0} ^ {k - 1} \beta_ {j} f _ {n + j} ^ {[ M ]}, \quad 0 \leqslant m \leqslant M - 1, \\ E: f _ {n + k} ^ {[ M ]} = f \left(t _ {n + k}, u _ {n + k} ^ {[ M ]}\right). \\ \end{array}
$$

原则上任一显式多步法和隐式多步法都可搭配成预-校算法及各种计算方案.

例1.1 Adams四阶四步预-校算法.取四阶四步Adams外插法为预估算法，四阶三步Adams内插法为校正算法，即得

$$
\begin{array}{l} P: u _ {n + 4} - u _ {n + 3} = \frac {h}{2 4} (5 5 f _ {n + 3} - 5 9 f _ {n + 2} + 3 7 f _ {n + 1} - 9 f _ {n}), \\ C: u _ {n + 4} - u _ {n + 3} = \frac {h}{2 4} (9 f _ {n + 4} + 1 9 f _ {n + 3} - 5 f _ {n + 2} + f _ {n + 1}). \\ \end{array}
$$

例1.2 Milne方法.以四阶四步法

$$
P: u _ {n + 4} - u _ {n} = \frac {4 h}{3} (2 f _ {n + 3} - f _ {n + 2} + 2 f _ {n + 1})
$$

为预估算法，四阶二步法

$$
C: u _ {n + 4} - u _ {n + 2} = \frac {h}{3} \left(f _ {n + 4} + 4 f _ {n + 3} + f _ {n + 2}\right)
$$

为校正算法, 得到由 $P$ 和 $C$ 组成的预-校方案 $PECE$ , 称为Milne算法, 计算公式为

$$
\begin{array}{l} P: \quad u _ {n + 4} ^ {[ 0 ]} - u _ {n} ^ {[ 1 ]} = \frac {4 h}{3} \left(2 f _ {n + 3} ^ {[ 1 ]} - f _ {n + 2} ^ {[ 1 ]} + 2 f _ {n + 1} ^ {[ 1 ]}\right), \\ E: \quad f _ {n + 4} ^ {[ 0 ]} = f \left(t _ {n + 4}, u _ {n + 4} ^ {[ 0 ]}\right), \\ C: \quad u _ {n + 4} ^ {[ 1 ]} - u _ {n + 2} ^ {[ 1 ]} = \frac {h}{3} \left(f _ {n + 4} ^ {[ 0 ]} + 4 f _ {n + 3} ^ {[ 1 ]} + f _ {n + 2} ^ {[ 1 ]}\right), \\ E: \quad f _ {n + 4} ^ {[ 1 ]} = f \left(t _ {n + 4}, u _ {n + 4} ^ {[ 1 ]}\right). \\ \end{array}
$$

# 1.2.4 多步法的计算问题

用 $k$ 步法计算时，需要知道 $k$ 个初值 $u_{0}, u_{1}, \dots, u_{k-1}$ ，其中 $u_{0} = u(t_{0})$ 是给定的初值，其余是附加初值。计算附加初值主要是用单步法，比如 Euler 法和 1.4 节将要介绍的 Runge-Kutta 法及其他单步法。为了保持多步法的精度，计算附加初值时要将 $t_{0}, t_{k-1}$ 之间的节点加密或采用和多步法有同样阶的 Runge-Kutta 法。

多步法的第二个问题是如何选择阶 $p$ (或者步数 $k$ ). 从收敛阶的观点, 自然希望把 $p$ 取大一些. 但是高阶收敛方法要求解的光滑性也高, 否则达不到高精度的目的. 从后面关于绝对稳定性的分析还知道, 高阶多步法的绝对稳定域也小, 所以 $p$ 的选取要考虑到解的光滑性和稳定性以及总的工作量.

多步法的第三个问题是步长 $h$ 的选取. 理论上似乎按照下节的误差估计式选定 $h$ 是合理的, 但那种估计往往偏大, 因此选定的 $h$ 可能过小, 既不必要也不经济. 实际用的步长 $h$ 不是一次取定, 而是根据精度要求, 由粗到细逐渐调整 (选步长), 当 $h$ 达到要求后就以此为步长计算. 在计算中还可以改变步长, 但计算过程变复杂了, 这是多步法的缺点. 与此相反, 单步法 (如 Euler 法及后面要介绍的 Runge-Kutta 法) 则适合用变步长计算.

# 1.2.5 习题

1. 用待定系数法求四阶三步方法类，确定四阶三步显式法  
2.满足条件 $\beta_{j} = 0,j = 0,1,2,\dots ,k - 1$ 的 $k$ 阶 $k$ 步法叫做Gear法，试对 $k = 1$ 2.3.4求Gear法的表达式

3. 用三阶 Adams 内插法及外插法分别解初值问题 $u' = -5u, u(0) = 1$ . 取步长 $h = 0.1$ , 0.05. 观察解在 $t = 1$ 处的误差，并与用 Euler 法计算的结果比较（参看 1.1 节习题 1）

# 1.3 相容性、稳定性和误差估计

本节讨论线性多步法的几个基本理论问题：相容性、稳定性和误差估计.

# 1.3.1 局部截断误差和相容性

考虑初值问题

$$
u ^ {\prime} = f (t, u), \quad t \in [ t _ {0}, T ] = [ 0, T ], \tag {1.70a}
$$

$$
u \left(t _ {0}\right) = u _ {0} \quad \text {题 向 美 门 的 长 时 间} \tag {1.70b}
$$

和逼近它的 $p$ 阶 $k$ 步法：

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u _ {n + j} = h \sum_ {j = 0} ^ {k} \beta_ {j} f _ {n + j}, \quad n = 0, 1 \dots . \tag {1.71}
$$

要想 (1.71) 的解 $u_{n}$ 逼近精确解 $u(t_{n})$ ，必需(1.71)在某种意义下逼近(1.70a).引进差分算子

$$
L [ u (t); h ] = \sum_ {j = 0} ^ {k} \left[ \alpha_ {j} u (t + j h) - h \beta_ {j} u ^ {\prime} (t + j h) \right]. \tag {1.72}
$$

设 $u(t)$ 是(1.70a)的具有 $p + 2$ 阶连续微商的解 $\left(u(t)\in C^{p + 2}\right)$ ，则由(1.62)和(1.63a)—(1.63b)，我们有

$$
\begin{array}{l} L [ u (t _ {n}); h ] = \sum_ {j = 0} ^ {k} \alpha_ {j} u (t _ {n} + j h) - h \sum_ {j = 0} ^ {k} \beta_ {j} u ^ {\prime} (t _ {n} + j h) \\ = c _ {p + 1} h ^ {p + 1} u ^ {(p + 1)} \left(t _ {n}\right) + O \left(h ^ {p + 2}\right) \tag {1.73} \\ \end{array}
$$

或

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u (t _ {n + j}) = h \sum_ {j = 0} ^ {k} \beta_ {j} f (t _ {n + j}, u (t _ {n + j})) + L [ u (t _ {n}); h ], \tag {1.74}
$$

$$
L \left[ u \left(t _ {n}\right); h \right] = c _ {p + 1} h ^ {p + 1} u ^ {(p + 1)} \left(t _ {n}\right) + O \left(h ^ {p + 2}\right), \tag {1.75}
$$

其中（参看(1.60))

$$
c _ {p + 1} = \frac {1}{(p + 1) !} \left(\alpha_ {1} + 2 ^ {p + 1} \alpha_ {2} + \dots + k ^ {p + 1} \alpha_ {k}\right) - \frac {1}{p !} \left(\beta_ {1} + 2 ^ {p} \beta_ {2} + \dots + k ^ {p} \beta_ {k}\right). \tag {1.76}
$$

像1.2节那样称 $L[u(t_n);h]$ 为局部截断误差，而称 $c_{p + 1}h^{p + 1}u^{(p + 1)}(t_n)$ 为局部截断误差的主项， $c_{p + 1}$ 为误差主项系数.在（1.74）中舍去 $L[u(t_n);h]$ ，并用 $u_{n + j}$ 代 $u(t_{n + j})$ 就导致多步法(1.71).我们关心的是误差 $e_n = u(t_n) - u_n$ ，称为整体误差.

现在考虑一般的 $k$ 步法 (1.71) (不必要求是 $p$ 阶方法). 为使 (1.71) 的解 $u_{n}$ 当 $h \to 0$ 时有可能收敛到 (1.70a) 的解 $u(t)$ , 自然要求

$$
\frac {1}{h} \left[ \sum_ {j = 0} ^ {k} \alpha_ {j} u (t _ {n + j}) - h \sum_ {j = 0} ^ {k} \beta_ {j} f (t _ {n + j}, u (t _ {n + j})) \right] - [ u ^ {\prime} (t _ {n}) - f (t _ {n}, u (t _ {n}) ] = o (1) (h \rightarrow 0), \tag {1.77}
$$

而 $u^{\prime}(t) = f(t,u(t))$ ，在(1.72）中令 $t = t_n$ ，则(1.77）可写成：

$$
L [ u (t _ {n}); h ] = o (h) (h \rightarrow 0). \tag {1.78}
$$

称多步法 (1.71) 相容, 如果对 (1.70a)(1.70b) 的任意光滑解 $u(t)$ , 关系 (1.77) 或 (1.78) 成立. 注意当 $u(t)$ 和 $f(t,u)$ 连续可微时, (1.77) 右端的 $o(1) = O(h)$ , 从而 $L[u(t);h] = O\left(h^{2}\right)$ , 所以多步法 (1.71) 至少是一阶的 (参看 (1.63a)(1.63b)). 这样可将相容性定义为

定义1.1 解初值问题(1.70a)(1.70b)的多步法(1.71)说是相容的（consistent)，如果它至少是一阶的。

引进多步法（1.71）的第一和第二特征多项式：

$$
\rho (\lambda) = \sum_ {j = 0} ^ {k} \alpha_ {j} \lambda^ {j}, \tag {1.79}
$$

$$
\sigma (\lambda) = \sum_ {j = 0} ^ {k} \beta_ {j} \lambda^ {j}. \tag {1.80}
$$

由(1.60)推出

定理1.1 为使 $k$ 步法（1.71）相容，必须且只需

$$
\rho (1) = 0, \quad \rho^ {\prime} (1) = \sigma (1). \tag {1.81}
$$

# 1.3.2 稳定性

用多步法计算时, 各种因素 (如初值 $u_{0}, u_{1}, \dots, u_{k-1}$ ) 是有误差的, 且这些误差将在计算中传递下去. 如果误差积累无限增长, 将会歪曲精确解, 这样的算法是不能用的, 为此我们对多步法提出稳定性要求.

定义1.2 称多步法(1.71)稳定(stable)，若存在常数 $C$ （不依赖 $h$ 和（1.71）的解）和 $h_0 > 0$ ，使 $\forall h\in (0,h_0)$ 和(1.71)的任何两个解 $\{u_n\}$ 和 $\{v_{n}\}$ （初值不同），恒有

$$
\max  _ {n h \leqslant T} | u _ {n} - v _ {n} | \leqslant C \max  _ {0 \leqslant j <   k} | u _ {j} - v _ {j} |. \tag {1.82}
$$

这等于说, 对一切充分小的 $h$ , 多步法的解连续依赖初值

定理1.2 设 $\rho(\lambda)$ 是形如(1.79)的第一特征多项式，则线性多步法(1.71)稳定的充要条件是 $\rho(\lambda)$ 满足根条件，即 $\rho(\lambda)$ 的所有根在单位圆内 $(|\lambda| \leqslant 1)$ ，且位于单位圆周上的根都是单根。

证明 必要性 将多步法用于方程 $u' = 0 (f = 0)$ . 此时 (1.71) 简化为

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u _ {n + j} = 0, \quad n = 0, 1, \dots ,
$$

其通解形如(1.19). 又 $\{v_{n} = 0\}$ 是上述方程的平凡解, 不等式 (1.82) 化为

$$
\max  _ {n h \leqslant T} | u _ {n} | \leqslant C \max  _ {0 \leqslant j <   k} | u _ {j} |, 0 <   h <   h _ {0},
$$

即 $\{u_n\}$ 关于 $n$ 和 $h$ ( $nh \leqslant T, 0 < h < h_0$ ) 一致有界. 而当 $h \to 0$ 时, $n$ 可趋于 $\infty$ , 由引理1.1的(i), $\rho(\lambda)$ 必满足根条件.

充分性 设 $\{u_n\}$ 和 $\{v_n\}$ 是 (1.71) 的任何两个解, 则

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u _ {n + j} = h \sum_ {j = 0} ^ {k} \beta_ {j} f (t _ {n + j}, u _ {n + j}),
$$

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} v _ {n + j} = h \sum_ {j = 0} ^ {k} \beta_ {j} f (t _ {n + j}, v _ {n + j}).
$$

令 $e_n = u_n - v_n$ ，则 $e_n$ 满足

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} e _ {n + j} = h b _ {n}, \quad n = 0, 1, \dots , \frac {T}{h}, \quad 0 <   h <   h _ {0}, \tag {1.83}
$$

其中

$$
b _ {n} = \sum_ {j = 0} ^ {k} \beta_ {j} [ f (t _ {n + j}, u _ {n + j}) - f (t _ {n + j}, v _ {n + j}) ]. \tag {1.84}
$$

设 $B = \max \left\{|\beta_0|, |\beta_1|, \dots, |\beta_k|\right\}$ , $f$ 关于 $u$ 满足 Lipschitz 条件:

$$
\left| f (t, u) - f (t, v) \right| \leqslant L | u - v |,
$$

则

$$
\left| b _ {n} \right| \leqslant B L \sum_ {j = 0} ^ {k} \left| e _ {n + j} \right|. \tag {1.85}
$$

引进向量 $E_{n} = (e_{n + k - 1},e_{n + k - 2},\dots ,e_{n})^{\mathrm{T}},B_{n} = (h\alpha_{k}^{-1}b_{n},0,\dots ,0)^{\mathrm{T}}(k$ 维）和矩阵

$$
C = \left[ \begin{array}{c c c c c} - \alpha_ {k} ^ {- 1} \alpha_ {k - 1} & - \alpha_ {k} ^ {- 1} \alpha_ {k - 2} & \dots & - \alpha_ {k} ^ {- 1} \alpha_ {1} & - \alpha_ {k} ^ {- 1} \alpha_ {0} \\ 1 & 0 & \dots & 0 & 0 \\ 0 & 1 & \dots & 0 & 0 \\ \vdots & \vdots & & \vdots & \vdots \\ 0 & 0 & \dots & 1 & 0 \end{array} \right]
$$

便可将 (1.83) 写成向量形式:

$$
\boldsymbol {E} _ {n + 1} = \boldsymbol {C E} _ {n} + \boldsymbol {B} _ {n},
$$

进而有

$$
\boldsymbol {E} _ {n} = C ^ {n} \boldsymbol {E} _ {0} + \sum_ {l = 0} ^ {n - 1} C ^ {l} \boldsymbol {B} _ {n - l - 1}, \quad n = 1, 2, \dots , \frac {T}{h}, \quad 0 <   h <   h _ {0}. \tag {1.86}
$$

今设 $\rho (\lambda)$ 满足根条件, 则由引理1.1的(i), 矩阵 $\{C^n\}$ 一致有界. 以 $\| E_n\|$ 表示向量的欧氏模, $\| C\|$ 表示相应的矩阵模, 则有常数 $M$ 使

$$
\| C ^ {n} \| \leqslant M, \quad n = 1, 2, \dots , \frac {T}{h}, \quad 0 <   h <   h _ {0}. \tag {1.87}
$$

又

$$
\left\| \boldsymbol {B} _ {n} \right\| \leqslant h \left| \alpha_ {k} ^ {- 1} \right| | b _ {n} | \leqslant B L \left| \alpha_ {k} ^ {- 1} \right| h \sum_ {j = 0} ^ {k} | e _ {n + j} |, \tag {1.88}
$$

于是由 (1.86) 和 (1.87) 得

$$
\| \boldsymbol {E} _ {n} \| \leqslant M \| \boldsymbol {E} _ {0} \| + M B L \left| \alpha_ {k} ^ {- 1} \right| h \sum_ {l = 0} ^ {n - 1} \sum_ {j = 0} ^ {k} | e _ {n + j - l - 1} |, \tag {1.89}
$$

而

$$
\begin{array}{l} \sum_ {l = 0} ^ {n - 1} \sum_ {j = 0} ^ {k} | e _ {n + j - l - 1} | = \sum_ {l = 0} ^ {n - 1} | e _ {n + k - l - 1} | + \sum_ {l = 0} ^ {n - 1} \sum_ {j = 0} ^ {k - 1} | e _ {n + j - l - 1} | \\ \leqslant \sum_ {l = 0} ^ {n - 1} \| E _ {n - l} \| + \sqrt {k} \sum_ {l = 0} ^ {n - 1} \| E _ {n - l - 1} \| \\ \leqslant \| \boldsymbol {E} _ {n} \| + (\sqrt {k} + 1) \sum_ {l = 0} ^ {n - 1} \| \boldsymbol {E} _ {n - l - 1} \| \\ \leqslant \| \boldsymbol {E} _ {n} \| + (\sqrt {k} + 1) \sum_ {j = 0} ^ {n - 1} \| \boldsymbol {E} _ {j} \|. \\ \end{array}
$$

故

$$
\| \boldsymbol {E} _ {n} \| \leqslant M \| \boldsymbol {E} _ {0} \| + M B L \left| \alpha_ {k} ^ {- 1} \right| h \| \boldsymbol {E} _ {n} \| + M B L \left| \alpha_ {k} ^ {- 1} \right| h (\sqrt {k} + 1) \sum_ {j = 0} ^ {n - 1} \| \boldsymbol {E} _ {j} \|.
$$

取 $h > 0$ 充分小, 例如 $h < h_0$ , 使

$$
M B L \left| \alpha_ {k} ^ {- 1} \right| h <   1.
$$

令

$$
\begin{array}{l} K _ {1} = \left(1 - M B L \left| \alpha_ {k} ^ {- 1} \right| h\right) ^ {- 1} M, \\ K _ {2} = \left(1 - M B L \left| \alpha_ {k} ^ {- 1} \right| h\right) ^ {- 1} M B L \left| \alpha_ {k} ^ {- 1} \right| (\sqrt {k} + 1), \\ \end{array}
$$

则

$$
\| \boldsymbol {E} _ {n} \| \leqslant K _ {1} \| \boldsymbol {E} _ {0} \| + K _ {2} h \sum_ {j = 0} ^ {n - 1} \| \boldsymbol {E} _ {j} \|. \tag {1.90}
$$

最后，利用 Gronwall 不等式 (1.26) 得到

$$
\left\| \boldsymbol {E} _ {n} \right\| \leqslant \mathrm {e} ^ {K _ {2} T} \left(K _ {1} + K _ {2} k h\right) \left\| \boldsymbol {E} _ {0} \right\|, \quad n = 1, 2, \dots , \frac {T}{h}, \quad 0 <   h <   h _ {0}. \tag {1.91}
$$

这证明了多步法稳定

对单步法 $(k = 1):\alpha_{1}u_{n + 1} + \alpha_{0}u_{n} = h(\beta_{1}f_{n + 1} + \beta_{0}f_{n}),\rho (\lambda) = \alpha_{1}\lambda +\alpha_{0},$ 唯一根 $\lambda_1 = -\frac{\alpha_0}{\alpha_1}$ 如果方法相容,则 $\alpha_{1} + \alpha_{0} = 0,\lambda_{1} = 1$ ，且是单根，故稳定.特别地,Euler法稳定.

对 Adams 外插法和内插法, 相应的 $\rho(\lambda) = \lambda^k - \lambda^{k-1} = (\lambda - 1)\lambda^{k-1}(k \geqslant 2)$ , 除单根 $\lambda = 1$ 在单位圆周上外, 其余的重根 $\lambda = 0$ 都在单位圆内部, 所以 Adams 法稳定. 若 $\rho(\lambda) = \lambda^k - \lambda^{k-2} = (\lambda^2 - 1)\lambda^{k-2}(k \geqslant 2)$ , 则称相应的显方法为 Nystrom 法, 相应的隐式方法为广义 Milne 法. 因为唯一可能的重根 $\lambda = 0$ 在单位圆内部, 而在单位圆周上的根 $\lambda = \pm 1$ 都是单根, 所以稳定. 特别地, Milne 方法 (1.65) 稳定. 不稳定的方法是不能用的.

例1.3 初值问题

$$
\begin{array}{l} u ^ {\prime} = 4 t u ^ {\frac {1}{2}}, \quad 0 \leqslant t \leqslant 2, \\ u (0) = 1 \\ \end{array}
$$

的精确解为 $u(t) = \left(1 + t^2\right)^2$ 考虑线性二步法：

$$
u _ {n + 2} - (1 + a) u _ {n + 1} + a u _ {n} = \frac {1}{2} h [ (3 - a) f _ {n + 1} - (1 + a) f _ {n} ],
$$

当 $a \neq -5$ 时是二阶方法, $a = -5$ 时是三阶方法. 第一特征多项式

$$
\rho (\lambda) = \lambda^ {2} - (1 + a) \lambda + a = (\lambda - 1) (\lambda - a).
$$

当 $a = 0$ 时稳定, $a = -5$ 时不稳定. 取步长 $h = 0.1$ , 初值 $u_0 = 1$ , 附加初值 $u_1 = \left(1 + h^2\right)^2 (h = 0.1)$ 是精确的. 用方案 (i) $a = 0$ 和 (ii) $a = -5$ 计算. 在刚开始的几步, 两种方案算出的结果都和精确解符合, 且 (ii) 比 (i) 更精确. 但再往后算, 方案 (i) 的结果仍

和精确解基本符合, 方案 (ii) 的误差则急剧增长, 完全歪曲了精确解, 具体数据如表 1.5.

表 1.5 具体数据  

<table><tr><td>t</td><td>精确解</td><td>(i) a = 0</td><td>(ii) a = -5</td></tr><tr><td>0.0</td><td>1.000 000 0</td><td>1.000 000 0</td><td>1.000 000 0</td></tr><tr><td>0.1</td><td>1.020 100 0</td><td>1.020 100 0</td><td>1.020 100 0</td></tr><tr><td>0.2</td><td>1.081 600 0</td><td>1.080 700 0</td><td>1.081 200 0</td></tr><tr><td>0.3</td><td>1.188 100 0</td><td>1.185 248 1</td><td>1.189 238 5</td></tr><tr><td>0.4</td><td>1.345 600 0</td><td>1.339 629 8</td><td>1.338 866 0</td></tr><tr><td>0.5</td><td>1.562 500 0</td><td>1.552 090 0</td><td>1.592 993 5</td></tr><tr><td>:</td><td>:</td><td>:</td><td>:</td></tr><tr><td>1.0</td><td>4.000 000 0</td><td>3.940 690 3</td><td>-68.639 804</td></tr><tr><td>1.1</td><td>4.884 100 0</td><td>4.808 219 7</td><td>367.263 92</td></tr><tr><td>:</td><td>:</td><td>:</td><td>:</td></tr><tr><td>2.0</td><td>25.000 000 0</td><td>24.632 457</td><td>-6.96 × 10^8</td></tr></table>

# 1.3.3 收敛性和误差估计

有了前述准备, 我们就可证明数值解的收敛性并估计整体误差. 设 $u(t)$ 为初值问题 (1.70a) (1.70b) 的解, $u_{n}$ 是线性多步法 (1.71) 的解, 则 $u(t_{n})$ 和 $u_{n}$ 分别满足

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u (t _ {n + j}) = h \sum_ {j = 0} ^ {k} \beta_ {j} f (t _ {n + j}, u (t _ {n + j})) + L [ u (t _ {n}); h ] \tag {1.92}
$$

和

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u _ {n + j} = h \sum_ {j = 0} ^ {k} \beta_ {j} f _ {n + j} \quad (f _ {n + j} = f (t _ {n + j}, u _ {n + j})). \tag {1.93}
$$

两式相减，则误差（整体误差） $e_n = u(t_n) - u_n$ 满足

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} e _ {n + j} = h b _ {n} + L [ u (t _ {n}); h ], \tag {1.94}
$$

$$
b _ {n} = \sum_ {j = 0} ^ {k} \beta_ {j} \left[ f \left(t _ {n + j}, u \left(t _ {n + j}\right)\right) - f \left(t _ {n + j}, u _ {n + j}\right) \right]. \tag {1.95}
$$

如前所设， $f$ 关于 $\mathcal{U}$ 满足Lipschitz条件，因此 $b_{n}$ 仍满足不等式(1.85).若 $k$ 步法(1.71）相容，例如，设（1.71）是 $\mathcal{P}$ 阶方法，则 $L[u(t_n);h]$ 有渐近表示(1.75)，于是可将(1.94）写成：

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} e _ {n + j} = h b _ {n} ^ {*}, \tag {1.96}
$$

其中

$$
b _ {n} ^ {*} = b _ {n} + c _ {p + 1} h ^ {p} u ^ {(p + 1)} \left(t _ {n}\right) + O \left(h ^ {p + 1}\right).
$$

令 $B_{n}^{*} = \left(h\alpha_{k}^{-1}b_{n}^{*},0,\dots ,0\right)^{\mathrm{T}}$ ，则由（1.85），有

$$
\begin{array}{l} \left\| \boldsymbol {B} _ {n} ^ {*} \right\| \leqslant \left| h \alpha_ {k} ^ {- 1} b _ {n} ^ {*} \right| \leqslant h \left| \alpha_ {k} ^ {- 1} \right| | b _ {n} | + M _ {p + 1} h ^ {p + 1} \\ \leqslant h B L \left| \alpha_ {k} ^ {- 1} \right| \sum_ {j = 0} ^ {k} | e _ {n + j} | + M _ {p + 1} h ^ {p + 1} \quad \left(M _ {p + 1} = \left| \alpha_ {k} ^ {- 1} \right| C _ {p + 1} \sup  _ {t} \left| u ^ {(p + 1)} (t) \right| + 1\right). \\ \end{array}
$$

比较 (1.96) 和 (1.83), 知向量 $\pmb{E}_{n} = (e_{n + k - 1}, e_{n + k - 2}, \dots, e_{n})^{\mathrm{T}}$ 仍可表示成 (1.86), 只需用 $\pmb{B}_{n}^{*}$ 代替那里的 $\pmb{B}_{n}$ .

若 $k$ 步法 (1.71) 稳定, 即 $\rho (\lambda)$ 满足根条件, 则 (1.87) 成立, 于是和 (1.89) 平行地有

$$
\| \boldsymbol {E} _ {n} \| \leqslant M \left(\left| \left| \boldsymbol {E} _ {0} \right| \right| + M _ {p + 1} T h ^ {p}\right) + M B L \left| \alpha_ {k} ^ {- 1} \right| h \sum_ {l = 0} ^ {n - 1} \sum_ {j = 0} ^ {k} | e _ {n + j - l - 1} |.
$$

取 $h > 0$ 充分小，使

$$
M B L \left| \alpha_ {k} ^ {- 1} \right| h <   1.
$$

则可得到与 (1.90) 平行的不等式:

$$
\left\| \boldsymbol {E} _ {n} \right\| \leqslant K _ {1} \left\| \boldsymbol {E} _ {0} + M _ {p + 1} T h ^ {p} \right\| + K _ {2} h \sum_ {j = 0} ^ {n - 1} \left\| \boldsymbol {E} _ {j} \right\|.
$$

最后由Gronwall不等式就得到误差估计：

$$
\left\| \boldsymbol {E} _ {n} \right\| \leqslant \mathrm {e} ^ {K _ {2} T} \left(\left(K _ {1} + K _ {2} k h\right) \left\| \boldsymbol {E} _ {0} \right\| + M _ {p + 1} T h ^ {p}\right). \tag {1.97}
$$

总之我们得

定理1.3 若解初值问题(1.70a)(1.70b)的多步法(1.71)相容而且稳定，则当 $h\to 0,t_n\to t$ 时数值解 $u_{n}\rightarrow u(t)$ ，其中 $\pmb{E}_0$ 是初始误差向量.若更设（1.71）是 $p$ 阶方法，则还有误差估计(1.97).

# 1.3.4 习题

1. 证明线性多步法

$$
u _ {n + 2} + (b - 1) u _ {n + 1} - b u _ {n} = \frac {1}{4} h [ (b + 3) f _ {n + 2} + (3 b + 1) f _ {n} ]
$$

当 $b \neq -1$ 时阶是2；当 $b = -1$ 时阶是3. 又 $b = -1$ 是不稳定的. 将 $b = -1$ 的方法用到 $u' = u, u(0) = 1$ ，解出相应的差分方程 $(u_0 = 1, u_1 = 1)$ ，说明方法发散

2. 确定 $\alpha$ 的变化域，使线性多步法

$$
u _ {n + 3} + \alpha (u _ {n + 2} - u _ {n + 1}) - u _ {n} = \frac {1}{2} (3 + \alpha) h (f _ {n + 2} + f _ {n + 1})
$$

是稳定的，并说明方法的阶不能大于2

# 1.4 单步法和 Runge-Kutta 法

Euler法是最简单的单步法. 单步法不需要附加初值, 所需的存储量小, 改变步长灵活, 但线性单步法的阶最多是2. 本节将介绍非线性(关于 $f$ ) 高阶单步法, 重点是Runge-Kutta法.

# 1.4.1 Taylor展开法

设初值问题

$$
\left\{ \begin{array}{l} u ^ {\prime} = f (t, u), \\ u \left(t _ {0}\right) = u _ {0} \end{array} \right.
$$

的解充分光滑. 将 $u(t)$ 在 $t_0$ 处用Taylor公式展开：

$$
u \left(t _ {1}\right) = u \left(t _ {0}\right) + h u ^ {\prime} \left(t _ {0}\right) + \frac {h ^ {2}}{2 !} u ^ {(2)} \left(t _ {0}\right) + \dots + \frac {h ^ {p}}{p !} u ^ {(p)} \left(t _ {0}\right) + O \left(h ^ {p + 1}\right), \tag {1.98}
$$

其中 $u(t_0) = u_0, u'(t_0) = f(t_0, u(t_0)) = f(t_0, u_0)$ ,

$$
\left\{ \begin{array}{l} u ^ {(2)} (t _ {0}) = \frac {\mathrm {d}}{\mathrm {d} t} f \Big | _ {t = t _ {0}} = [ f _ {t} ^ {\prime} + f _ {u} ^ {\prime} \cdot u ^ {\prime} (t) ] _ {t = t _ {0}} \\ \quad = f _ {t} ^ {\prime} (t _ {0}, u _ {0}) + f (t _ {0}, u _ {0}) f _ {u} ^ {\prime} (t _ {0}, u _ {0}), \\ u ^ {(3)} (t _ {0}) = \frac {\mathrm {d}}{\mathrm {d} t} \left[ \frac {\mathrm {d}}{\mathrm {d} t} f \right] _ {t = t _ {0}} \\ \quad = f _ {t t} ^ {\prime \prime} (t _ {0}, u _ {0}) + 2 f (t _ {0}, u _ {0}) f _ {t u} ^ {\prime \prime} (t _ {0}, u _ {0}) + \\ \quad f ^ {2} (t _ {0}, u _ {0}) f _ {u u} ^ {\prime \prime} (t _ {0}, u _ {0}) + f _ {t} ^ {\prime} (t _ {0}, u _ {0}) f _ {u} ^ {\prime} (t _ {0}, u _ {0}) + \\ \quad f (t _ {0}, u _ {0}) f _ {u} ^ {\prime 2} (t _ {0}, u _ {0}), \\ \quad \dots \dots \dots \dots . \end{array} \right. \tag {1.99}
$$

令

$$
\varphi (t, u (t), h) = \sum_ {j = 1} ^ {p} \frac {h ^ {j - 1}}{j !} \frac {\mathrm {d} ^ {j - 1}}{\mathrm {d} t ^ {j - 1}} f (t, u (t)), \tag {1.100}
$$

则可将 (1.98) 改写为

$$
u \left(t _ {0} + h\right) - u \left(t _ {0}\right) = h \varphi \left(t _ {0}, h \left(t _ {0}\right), h\right) + O \left(h ^ {p + 1}\right). \tag {1.101}
$$

舍去余项 $O\left(h^{p + 1}\right)$ ，则得

$$
u _ {1} - u _ {0} = h \varphi (t _ {0}, u _ {0}, h).
$$

一般说来，若已知 $u_{n}$ ，则

$$
u _ {n + 1} - u _ {n} = h \varphi (t _ {n}, u _ {n}, h), \quad n = 0, 1, \dots . \tag {1.102}
$$

这是一个单步法, 局部截断误差为 $O(h^{p + 1})$ . 由 (1.99) 和 (1.100) 可知 $\varphi$ 关于 $f$ 非线性. 当 $p = 1$ 时它是 Euler 折线法. 由于计算 $\varphi(t_{n}, u_{n}, h)$ 的工作量太大, 一般不用 Taylor 展开法作数值计算, 但可用它计算附加初值.

# 1.4.2 单步法的稳定性和收敛性

将初值问题写成积分形式：

$$
u (t + h) - u (t) = \int_ {t} ^ {t + h} f (\tau , u (\tau)) \mathrm {d} \tau , \quad u \left(t _ {0}\right) = u _ {0}. \tag {1.103}
$$

如果有某一确定的函数 $\varphi (t,u,h)$ （通过某种离散化)，使初值问题的任一解 $u(t)$ 满足

$$
u (t + h) - u (t) = h \varphi (t, u (t), h) + O \left(h ^ {p + 1}\right), \tag {1.104}
$$

其中 $p \geqslant 1$ 是使 (1.104) 成立的最大整数, 则称算法

$$
u _ {n + 1} = u _ {n} + h \varphi (t _ {n}, u _ {n}, h), \quad n = 0, 1, \dots , \frac {T}{h} \tag {1.105}
$$

为 $p$ 阶单步法

Taylor展开法是 $p$ 阶单步法， $\varphi (t,u,h)$ 由(1.100)定义.Euler法是一阶单步法，相应的 $\varphi = f(t,u(t))$

注意

$$
\varphi (t, u (t), h) = \frac {u (t + h) - u (t)}{h} + O \left(h ^ {p}\right) = \frac {1}{h} \int_ {t} ^ {t + h} f (\tau , u (\tau)) \mathrm {d} \tau + O \left(h ^ {p}\right),
$$

故

$$
\lim  _ {h \rightarrow 0} \varphi (t, u (t), h) = f (t, u (t)).
$$

定义 $\varphi (t,u(t),0) = f(t,u(t))$ ，则知 $\varphi (t,u(t),h)$ 于 $h = 0$ 连续.反之，若 $\varphi (t,u(t),h)$ 于 $h = 0$ 连续且单步法(1.105）的局部截断误差 $R_{h} = o(h)(h\to 0)$ ，则由

$$
\varphi (t, u (t), h) = \frac {1}{h} \int_ {t} ^ {t + h} f (\tau , u (\tau)) \mathrm {d} \tau + o (1),
$$

令 $h\to 0$ ，知 $\varphi (t,u(t),0) = f(t,u(t))$ .所以可将单步法的相容性定义为

定义1.3 称单步法(1.105)相容，如果 $\varphi (t,u(t),h)$ 于 $h = 0$ 连续，且

$$
\varphi (t, u, 0) = f (t, u).
$$

定理1.4 设 $\varphi(t, u, h) (t_0 \leqslant t \leqslant T, 0 \leqslant h \leqslant h_0, u \in (-\infty, \infty))$ 关于 $u$ 满足 Lipschitz 条件，则单步法(1.105)稳定。

实际上, 设 $v_{n}$ 是 (1.105) 的以 $v_{0}$ 为初值的解, 则

$$
v _ {n + 1} = v _ {n} + h \varphi (t _ {n}, v _ {n}, h).
$$

与 (1.105) 相减, 知 $e_n = u_n - v_n$ 满足

$$
e _ {n + 1} = e _ {n} + h [ \varphi (t _ {n}, u _ {n}, h) - \varphi (t _ {n}, v _ {n}, h) ],
$$

于是

$$
\begin{array}{l} \left| e _ {n + 1} \right| \leqslant \left| e _ {n} \right| + h L \left| e _ {n} \right| = (1 + h L) \left| e _ {n} \right| \\ \leqslant \dots \leqslant (1 + h L) ^ {n + 1} | e _ {0} | \\ \leqslant \mathrm {e} ^ {L (n + 1) h} | e _ {0} | \leqslant \mathrm {e} ^ {L (T - t _ {0})} | e _ {0} |, \quad (n + 1) h \leqslant T - t _ {0}, \\ \end{array}
$$

所以（1.105）稳定

定理1.5 设 $\varphi(t, u, h)$ 满足定理1.4的条件，又单步法(1.105）相容，则当 $h \to 0$ 时，它的数值解 $u_n \to u(t)$ ，只要 $t_0 + nh \to t$ ，初值 $u_0 \to u(t_0)$ 。若更设(1.105）是 $p$ 阶单步法，则还有敛速估计(1.106)。

证明 不妨设 (1.105) 是 $p$ 阶单步法, 则初值问题的解 $u(t)$ 满足

$$
u \left(t _ {n + 1}\right) = u \left(t _ {n}\right) + h \varphi \left(t _ {n}, u \left(t _ {n}\right), h\right) + O \left(h ^ {p + 1}\right).
$$

与(1.105)相减, 并令 $e_n = u(t_n) - u_n$ ，则

$$
e _ {n + 1} = e _ {n} + h \left[ \varphi \left(t _ {n}, u \left(t _ {n}\right), h\right) - \varphi \left(t _ {n}, u _ {n}, h\right) \right] + O \left(h ^ {p + 1}\right).
$$

因 $\varphi$ 关于 $u$ 满足Lipschitz条件，所以

$$
\begin{array}{l} \left| e _ {n + 1} \right| \leqslant \left| e _ {n} \right| + h L \left| e _ {n} \right| + c h ^ {p + 1}, \\ \left| e _ {n + 1} \right| - \left| e _ {n} \right| \leqslant h L \left| e _ {n} \right| + c h ^ {p + 1}. \\ \end{array}
$$

两端关于 $n$ 求和，得

$$
\left| e _ {n} \right| \leqslant \left(\left| e _ {0} \right| + c T h ^ {p}\right) + h L \sum_ {k = 0} ^ {n - 1} \left| e _ {k} \right|.
$$

再利用Gronwall不等式(1.26)，就得到估计：

$$
\left| e _ {n} \right| \leqslant \mathrm {e} ^ {L T} \left[ (1 + L h) \left| e _ {0} \right| + c T h ^ {p} \right]. \tag {1.106}
$$

特别地, 若取 $u_0 = u(t_0)$ , 则 $e_0 = 0$ , 误差 $e_n$ 的阶为 $O(h^p)$ .

# 1.4.3 Runge-Kutta法

Taylor展开法，用 $f$ 在同一点 $(t_n,u_n)$ 的高阶导数表示 $\varphi (t_{n},u_{n},h)$ ，这不便于数值计算.Runge-Kutta法是用 $f$ 在一些点上的值表示 $\varphi (t_n,u_n,h)$ ，使单步法局部截断误差的阶和Taylor展开法相等.我们先在区间 $[t,t + h]$ 上讨论.将初值问题写成积分形式：

$$
u (t + h) = u (t) + \int_ {t} ^ {t + h} f (\tau , u (\tau)) \mathrm {d} \tau . \tag {1.107}
$$

在 $[t, t + h]$ 取 $m$ 个点 $t_1 = t \leqslant t_2 \leqslant t_3 \leqslant \dots \leqslant t_m \leqslant t + h$ . 若知道 $k_i = f(t_i, u(t_i)), i = 1, 2, \dots, m$ , 则可用它们的一次组合去近似 $f$ :

$$
\sum_ {i = 1} ^ {m} c _ {i} k _ {i} \approx f. \tag {1.108}
$$

问题是如何计算 $k_{i}$ (因 $u(t_{i})$ 未知). 一个直观的想法是: 设已知 $(t_{1}, k_{1}) = (t_{1}, f(t_{1}, u(t_{1})))$ , 由 Euler 法 $u(t_{2}) \approx u(t_{1}) + (t_{2} - t_{1})f(t_{1}, u(t_{1})) = u(t_{1}) + (t_{2} - t_{1})k_{1}$ , 于是

$$
k _ {2} \approx f (t _ {2}, u (t _ {1}) + (t _ {2} - t _ {1}) k _ {1}).
$$

再利用Euler法又可以由 $(t_2,k_2)$ 算出

$$
k _ {3} \approx f \left(t _ {3}, u \left(t _ {1}\right) + \left(t _ {2} - t _ {1}\right) k _ {1} + \left(t _ {3} - t _ {2}\right) k _ {2}\right).
$$

如此可继续下去. 要求节点 $\{t_i\}$ 和系数 $\{c_i\}$ 适当选取, 使近似式 (1.108) 有尽可能高的逼近阶. 这为下面构造 Runge-Kutta 法提供某些启示.

为便于推导，我们先引进若干记号. 首先令

$$
t _ {i} = t + a _ {i} h = t _ {1} + a _ {i} h, \quad i = 2, 3, \dots , m,
$$

其中 $a_{i}$ 与 $h$ 无关. 再引进下三角形系数阵:

$$
\begin{array}{c c c c} b _ {2 1} & & \\ b _ {3 1} & b _ {3 2} & & \\ \vdots & \vdots & \ddots & \\ b _ {m 1} & b _ {m 2} & \dots & b _ {m, m - 1}, \end{array}
$$

其中 $b_{ij}$ 与 $h$ 无关，

$$
\sum_ {j = 1} ^ {i - 1} b _ {i j} = a _ {i}, \quad i = 2, 3, \dots , m. \tag {1.109}
$$

又 $c_{i}\geqslant 0$

$$
\sum_ {i = 1} ^ {m} c _ {i} = 1. \tag {1.110}
$$

假设三组系数 $\{a_i\}, \{b_{ij}\}$ 和 $\{c_i\}$ 已给定, 则 Runge-Kutta 法计算过程如下:

$$
u _ {n + 1} = u _ {n} + h \varphi (t _ {n}, u _ {n}, h), \quad n = 0, 1, \dots , \tag {1.111}
$$

其中

$$
\varphi (t, u (t), h) = \sum_ {i = 1} ^ {m} c _ {i} k _ {i}, \tag {1.112}
$$

$$
\left\{ \begin{array}{l} k _ {1} = f (t, u), \\ k _ {2} = f (t + h a _ {2}, u (t) + h b _ {2 1} k _ {1}), \quad b _ {2 1} = a _ {2}, \\ k _ {3} = f (t + h a _ {3}, u (t) + h (b _ {3 1} k _ {1} + b _ {3 2} k _ {2})), \quad b _ {3 1} + b _ {3 2} = a _ {3}, \\ \dots \dots \dots \dots . \\ k _ {m} = f \left(t + h a _ {m}, u (t) + h \sum_ {j = 1} ^ {m - 1} b _ {m j} k _ {j}\right), \quad \sum_ {j = 1} ^ {m - 1} b _ {m j} = a _ {m}. \end{array} \right. \tag {1.113}
$$

系数 $\{a_i\}, \{b_{ij}\}$ 和 $\{c_i\}$ 按如下原则确定: 将 $k_i$ 关于 $h$ 展开, 以之代到 (1.112), 使 $l$ 次幂 $h^l (l = 0,1,\dots ,p - 1)$ 的系数和 (1.100) 同次幂的系数相等, 如此得到的算法 (1.111) 称为 $m$ 级 $p$ 阶 Runge-Kutta 法.

现在推导一些常用的计算方案. 将 $u(t + h)$ 展开到 $h$ 的三次幂:

$$
u (t + h) = u (t) + \sum_ {l = 1} ^ {3} \frac {h ^ {l}}{l !} u ^ {(l)} (t) + O \left(h ^ {4}\right) = u (t) + h \varphi_ {T} (t, u, h), \tag {1.114}
$$

其中

$$
\left\{ \begin{array}{l} \varphi_ {T} (t, u, h) = f + \frac {1}{2} h F + \frac {1}{6} h ^ {2} \left(F f _ {u} ^ {\prime} + G\right) + O \left(h ^ {3}\right), \\ F = f _ {t} ^ {\prime} + f f _ {u} ^ {\prime}, \\ G = f _ {t t} ^ {\prime \prime} + 2 f f _ {t u} ^ {\prime \prime} + f ^ {2} f _ {u u} ^ {\prime \prime}. \end{array} \right. \tag {1.115}
$$

其次，由二元Taylor展开式

$$
\begin{array}{l} k _ {1} = f (t, u) = f, \\ k _ {2} = f \left(t + h a _ {2}, u + h a _ {2} k _ {1}\right) \\ = f + h a _ {2} \left(f _ {t} ^ {\prime} + k _ {1} f _ {u} ^ {\prime}\right) + \frac {1}{2} h ^ {2} a _ {2} ^ {2} \left(f _ {t t} ^ {\prime \prime} + 2 k _ {1} f _ {t u} ^ {\prime \prime} + k _ {1} ^ {2} f _ {u u} ^ {\prime \prime}\right) + O \left(h ^ {3}\right) \\ \end{array}
$$

$$
= f + h a _ {2} F + \frac {1}{2} h ^ {2} a _ {2} ^ {2} G + O \left(h ^ {3}\right).
$$

同样地,

$$
k _ {3} = f + h a _ {3} F + h ^ {2} \left(a _ {2} b _ {3 2} f _ {u} ^ {\prime} F + \frac {1}{2} a _ {3} ^ {2} G\right) + O \left(h ^ {3}\right).
$$

于是

$$
\begin{array}{l} \varphi (t, u, h) = \left(c _ {1} + c _ {2} + c _ {3}\right) f + h \left(a _ {2} c _ {2} + a _ {3} c _ {3}\right) F + \\ \frac {1}{2} h ^ {2} \left[ 2 a _ {2} b _ {3 2} c _ {3} f _ {u} ^ {\prime} F + \left(a _ {2} ^ {2} c _ {2} + a _ {3} ^ {2} c _ {3}\right) G \right] + O \left(h ^ {3}\right). \tag {1.116} \\ \end{array}
$$

比较 $\varphi (t,u,h)$ 和 $\varphi_T(t,u,h)$ 的同次幂系数，可得以下具体方案

1. $m = 1$ 比较 $h$ 的零次幂，知

$$
\varphi (t, u, h) = f,
$$

算法 (1.111) 是 Euler 法

2. $m = 2$ 此时 $c_{3} = 0$

$$
\varphi (t, u, h) = \left(c _ {1} + c _ {2}\right) f + h a _ {2} c _ {2} F + \frac {1}{2} h ^ {2} a _ {2} ^ {2} c _ {2} G + O \left(h ^ {3}\right).
$$

与 $\varphi_T(t,u,h)$ 比较 $1,h$ 的系数，则

$$
c _ {1} + c _ {2} = 1, \quad a _ {2} c _ {2} = \frac {1}{2}.
$$

它有无穷多组解，从而有无穷多个二级二阶算法. 两个常见的方法是

(1) $c_{1} = 0, c_{2} = 1, a_{2} = \frac{1}{2}$ . 此时

$$
u _ {n + 1} = u _ {n} + h f \left(t _ {n} + \frac {1}{2} h, u _ {n} + \frac {1}{2} h f _ {n}\right).
$$

称为中点法, 这是一种修正的 Euler 法.

(2) $c_{1} = c_{2} = \frac{1}{2}, a_{2} = 1$ , 此时

$$
u _ {n + 1} = u _ {n} + \frac {1}{2} h (f (t _ {n}, u _ {n}) + f (t _ {n + 1}, u _ {n} + h f _ {n})).
$$

这是改进的Euler法

3. $m = 3$ 比较(1.115)和(1.116)，令 $1,h,h^2$ 的系数相等，并注意 $F,G$ 的任意性，得

$$
\begin{array}{l} c _ {1} + c _ {2} + c _ {3} = 1, \quad a _ {2} c _ {2} + a _ {3} c _ {3} = \frac {1}{2}, \\ a _ {2} ^ {2} c _ {2} + a _ {3} ^ {2} c _ {3} = \frac {1}{3}, \quad a _ {2} b _ {3 2} c _ {3} = \frac {1}{6}. \\ \end{array}
$$

由这4个方程不能完全确定6个系数，因此这是含两个参数的三级三阶方法.常见的

方案有

(1) Heun 三阶方法. 此时

$$
c _ {1} = \frac {1}{4}, \quad c _ {2} = 0, \quad c _ {3} = \frac {3}{4},
$$

$$
a _ {2} = \frac {1}{3}, \quad a _ {3} = \frac {2}{3}, \quad b _ {3 2} = \frac {2}{3}.
$$

算法为

$$
\left\{ \begin{array}{l} u _ {n + 1} = u _ {n} + \frac {h}{4} \left(k _ {1} + 3 k _ {3}\right), \\ k _ {1} = f \left(t _ {n}, u _ {n}\right), \\ k _ {2} = f \left(t _ {n} + \frac {1}{3} h, u _ {n} + \frac {1}{3} h k _ {1}\right), \\ k _ {3} = f \left(t _ {n} + \frac {2}{3} h, u _ {n} + \frac {2}{3} h k _ {2}\right). \end{array} \right. \tag {1.117}
$$

(2)Kutta三阶方法.此时

$$
c _ {1} = \frac {1}{6}, \quad c _ {2} = \frac {2}{3}, \quad c _ {3} = \frac {1}{6},
$$

$$
a _ {2} = \frac {1}{2}, \quad a _ {3} = 1, \quad b _ {3 2} = 2.
$$

算法为

$$
\left\{ \begin{array}{l} u _ {n + 1} = u _ {n} + \frac {h}{6} \left(k _ {1} + 4 k _ {2} + k _ {3}\right), \\ k _ {1} = f \left(t _ {n}, u _ {n}\right), \\ k _ {2} = f \left(t _ {n} + \frac {1}{2} h, u _ {n} + \frac {1}{2} h k _ {1}\right), \\ k _ {3} = f \left(t _ {n} + h, u _ {n} - h k _ {1} + 2 h k _ {2}\right). \end{array} \right. \tag {1.118}
$$

当 $f$ 与 $u$ 无关时，这就是Simpson公式

4. $m = 4$ 将(1.115)，(1.116)展开到 $h^3$ ，比较 $h^i (i = 0,1,2,3)$ 的系数，则得含13个待定系数的11个方程，由此得到含两个参数的四级四阶Runge-Kutta方法类，其中最常用的有以下两个算法：

四阶 Runge-Kutta 法:

$$
\left\{ \begin{array}{l} u _ {n + 1} = u _ {n} + \frac {h}{6} \left(k _ {1} + 2 k _ {2} + 2 k _ {3} + k _ {4}\right), \\ k _ {1} = f \left(t _ {n}, u _ {n}\right), \\ k _ {2} = f \left(t _ {n} + \frac {1}{2} h, u _ {n} + \frac {1}{2} h k _ {1}\right), \\ k _ {3} = f \left(t _ {n} + \frac {1}{2} h, u _ {n} + \frac {1}{2} h k _ {2}\right), \\ k _ {4} = f \left(t _ {n} + h, u _ {n} + h k _ {3}\right) \end{array} \right. \tag {1.119}
$$

和

$$
\left\{ \begin{array}{l} u _ {n + 1} = u _ {n} + \frac {h}{8} \left(k _ {1} + 3 k _ {2} + 3 k _ {3} + k _ {4}\right), \\ k _ {1} = f \left(t _ {n}, u _ {n}\right), \\ k _ {2} = f \left(t _ {n} + \frac {1}{3} h, u _ {n} + \frac {1}{3} h k _ {1}\right), \\ k _ {3} = f \left(t _ {n} + \frac {2}{3} h, u _ {n} - \frac {1}{3} h k _ {1} + h k _ {2}\right), \\ k _ {4} = f \left(t _ {n} + h, u _ {n} + h k _ {1} - h k _ {2} + h k _ {3}\right). \end{array} \right. \tag {1.120}
$$

(1.119) 是最常用的 Runge-Kutta 法, 通常也称为经典的四阶 Runge-Kutta 法.

仿此人们还可造出更多的算法. 设 $f(t,u)$ $(t_0\leqslant t\leqslant T,u\in (-\infty , + \infty))$ 连续, 且关于 $u$ 满足Lipschitz条件. 由 $\varphi (t,u(t),h)$ 的表达式(1.112)(1.113）（注意 $c_{1} + c_{2} + \dots +c_{m} =$ 1)知 $\varphi (t,u,0) = f(t,u)$ ，即Runge-Kutta法相容，所以定理1.4和定理1.5对Runge-Kutta法恒成立，特别地， $p$ 阶Runge-Kutta法整体误差的阶为 $O(h^{p})$

例1.4 用四级四阶Runge-Kutta法计算初值问题：

$$
u ^ {\prime} = 4 t u ^ {\frac {1}{2}}, \quad 0 \leqslant t \leqslant 2,
$$

$$
u (0) = 1.
$$

取 $h = 0.1,0.5,1$ 精确解为

$$
u (t) = \left(1 + t ^ {2}\right) ^ {2}.
$$

计算结果如表1.6.与例1.3的表1.5比较，可见Runge-Kutta法的稳定性较线性三阶二步法优越

表 1.6 四阶 Runge-Kutta 法的计算结果  

<table><tr><td>t</td><td>精确解</td><td>h=0.1</td><td>h=0.5</td><td>h=1</td></tr><tr><td>0.0</td><td>1.000 000</td><td>1.000 000</td><td>1.000 000</td><td>1.000 000</td></tr><tr><td>0.3</td><td>1.081 600</td><td>1.081 599</td><td>-</td><td>-</td></tr><tr><td>0.5</td><td>1.562 500</td><td>1.562 497</td><td>1.561 106</td><td>-</td></tr><tr><td>0.8</td><td>2.689 600</td><td>2.689 592</td><td>-</td><td>-</td></tr><tr><td>1.0</td><td>4.000 000</td><td>3.999 985</td><td>3.993 247</td><td>3.913 900</td></tr><tr><td>1.3</td><td>5.953 600</td><td>5.953 576</td><td>-</td><td>-</td></tr><tr><td>1.5</td><td>10.562 500</td><td>10.562 455</td><td>10.542 656</td><td>-</td></tr><tr><td>2.0</td><td>25.000 000</td><td>24.999 904</td><td>24.957 954</td><td>24.530 977</td></tr></table>

# 1.4.4 习题

就 $c_{2} = c_{3}$ 和 $a_2 = a_3$ 导出三阶Runge-Kutta法

# *1.5 绝对稳定性和绝对稳定域

# 1.5.1 绝对稳定性

无论从理论还是应用方面看，单步法和多步法都必须是稳定的。但这种稳定有两个限制，一是要求 $h \in (0, h_0)$ 充分小，而实际用的 $h$ 是固定的；二是只允许初值有误差，往后各步计算都精确，而实际计算时每步都可能有舍入误差。为了控制这种误差的增长，需对多步法提出进一步要求，即绝对稳定性。

对一般非线性常微分方程组：

$$
\boldsymbol {u} ^ {\prime} = \boldsymbol {f} (t, \boldsymbol {u})
$$

讨论绝对稳定性是困难的，通常是考虑它在解 $\pmb{u}$ 临近的线性化方程：

$$
\left(\boldsymbol {u} - \overline {{\boldsymbol {u}}}\right) ^ {\prime} = \frac {\partial f (t , \overline {{\boldsymbol {u}}})}{\partial \boldsymbol {u}} \left(\boldsymbol {u} - \overline {{\boldsymbol {u}}}\right).
$$

为此我们讨论线性常微分方程组：

$$
\boldsymbol {u} ^ {\prime} = \boldsymbol {A} \boldsymbol {u}, \tag {1.121}
$$

其中 $\mathbf{A}$ 是常数矩阵. 设 $\mathbf{A}$ 可对角化:

$$
\boldsymbol {T} ^ {- 1} \boldsymbol {A} \boldsymbol {T} = \operatorname {d i a g} (\lambda_ {1}, \lambda_ {2}, \dots , \lambda_ {m}),
$$

$\lambda_{j}$ 是实或复特征值，作变换 $\pmb {u} = \pmb{T}\pmb{v}$ ，则(1.121）化为

$$
v _ {i} ^ {\prime} = \lambda_ {i} v _ {i}, \quad i = 1, 2, \dots , m.
$$

因此作为模型，我们考虑方程：

$$
u ^ {\prime} = \mu u \tag {1.122}
$$

的绝对稳定性, 其中 $\mu$ 是实数或复数

求解 (1.122) 的线性多步法为

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u _ {n + j} = \mu h \sum_ {j = 0} ^ {k} \beta_ {j} u _ {n + j}. \tag {1.123}
$$

令 $\bar{h} = \mu h$ ，则

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} u _ {n + j} - \bar {h} \sum_ {j = 0} ^ {k} \beta_ {j} u _ {n + j} = 0. \tag {1.124}
$$

实际计算时每步都可能有舍入误差，得到的是 $u_{n}$ 的近似 $\bar{u}_n$ ，满足

$$
\sum_ {j = 0} ^ {k} \left(\alpha_ {j} - \bar {h} \beta_ {j}\right) \bar {u} _ {n + j} = \eta_ {n}, \tag {1.125}
$$

其中 $\eta_{n}$ 为局部舍入误差.假设 $|\eta_n|\leqslant M.$ 令 $\bar{e}_n = \bar{u}_n - u_n$ ，显然 $\bar{e}_n$ 满足

$$
\sum_ {j = 0} ^ {k} \left(\alpha_ {j} - \bar {h} \beta_ {j}\right) \bar {e} _ {n + j} = \eta_ {n}. \tag {1.126}
$$

如前引进 $k$ 维列向量 $\pmb{E}_{n} = (\bar{e}_{n + k - 1},\bar{e}_{n + k - 2},\dots ,\bar{e}_{n})^{\mathrm{T}}$ 和 $\pmb {\eta}_n = (\alpha_k^{-1}\eta_n,0,\dots ,0)^{\mathrm{T}}$ 以及 $k\times k$ 矩阵

$$
\overline {{C}} = \left[ \begin{array}{c c c c c} - a _ {k} ^ {- 1} a _ {k - 1} & - a _ {k} ^ {- 1} a _ {k - 2} & \dots & - a _ {k} ^ {- 1} a _ {1} & - a _ {k} ^ {- 1} a _ {0} \\ 1 & 0 & \dots & 0 & 0 \\ 0 & 1 & \dots & 0 & 0 \\ \vdots & \vdots & & \vdots & \vdots \\ 0 & 0 & \dots & 1 & 0 \end{array} \right], \tag {1.127}
$$

其中 $a_{j} = \alpha_{j} - \bar{h}\beta_{j}, j = k, k-1,\dots, 1,0$ . 假定 $h$ 充分小使 $a_{k} \neq 0$ , 则可将 (1.126) 写成向量形式:

$$
\boldsymbol {E} _ {n} = \overline {{C}} \boldsymbol {E} _ {n - 1} + \boldsymbol {\eta} _ {n}.
$$

逐次递推，得

$$
\boldsymbol {E} _ {n} = \overline {{\boldsymbol {C}}} ^ {n} \boldsymbol {E} _ {0} + \left(\overline {{\boldsymbol {C}}} ^ {n - 1} \eta_ {0} + \dots + \overline {{\boldsymbol {C}}} ^ {n - i - 1} \eta_ {i} + \dots + \eta_ {n - 1}\right) \tag {1.128}
$$

(参看 (1.86)). 为使误差 $E_{n}$ 在逐步计算中减小, 应要求

$$
\lim  _ {n \rightarrow \infty} \overline {{C}} ^ {n} = O. \tag {1.129}
$$

令

$$
\rho (\lambda) = \sum_ {j = 0} ^ {k} \alpha_ {j} \lambda^ {j}, \quad \sigma (\lambda) = \sum_ {j = 0} ^ {k} \beta_ {j} \lambda^ {j}, \tag {1.130}
$$

则矩阵 $\overline{C}$ 的特征方程为

$$
\rho (\lambda) - \bar {h} \sigma (\lambda) = 0. \tag {1.131}
$$

由引理1.1的(iii)，可知(1.129)成立的充要条件是(1.131)的根都在单位圆内

定义1.4 称线性 $k$ 步法关于 $\bar{h} = \mu h$ 绝对稳定（absolutely stable），如果特征方程(1.131)的根都在单位圆内 $(|\lambda| < 1)$ . 若存在复数域 $D_A$ ，使多步法对任意的 $\bar{h} \in D_A$ 都绝对稳定，则称 $D_A$ 为绝对稳定域.

显然绝对稳定域越大，方法的适应性越强，因此也更优越

可以证明，当矩阵 $\overline{C}$ 的特征值都在单位圆内时（ $\overline{C}$ 的谱半径小于1)，存在一种范数

$\| \cdot \|$ ，使 $\| \overline{C}\| < 1.$ 由(1.128)，直接得

$$
\| E _ {n} \| \leqslant \| \overline {{C}} \| ^ {n} + M (1 - \| \overline {{C}} \|) ^ {- 1}, \quad M \geqslant | \eta_ {n} |, \quad n = 1, 2, \dots ,
$$

所以舍入误差的影响是可控的 (连续依赖 $M$ ).

对一般的非线性问题， $f_{u}^{\prime}$ 不是常数，此时可取 $\mu$ 是 $f_{u}^{\prime}$ 的界，或是 $f_{u}^{\prime}$ 的一个或多个代表值， $h$ 为最大容许步长.

# 1.5.2 绝对稳定域

设 $D_A$ 是多步法的绝对稳定域, 则特征方程 (1.131) 确定一由单位圆 $|\lambda| < 1$ 到复平面 $\bar{h} \in Z$ 的解析变换:

$$
\bar {h} = \frac {\rho (\lambda)}{\sigma (\lambda)}. \tag {1.132}
$$

命题1.1 一个相容、稳定的多步法若绝对稳定，则 $\mu$ 的实部 $\operatorname{Re}\mu < 0$ ，从而绝对稳定域 $D_A \subseteq Z_{-}$ (左半平面).

证明 因多步法相容且稳定, 故 1 必为 $\rho(\lambda)$ 的单根, 即 $\rho(1) = 0$ , 但 $\sigma(1) = \rho'(1) \neq 0$ . 设 $k$ 次方程

$$
\rho (\lambda) - \bar {h} \sigma (\lambda) = 0
$$

的根为 $\lambda_1(\bar{h}),\lambda_2(\bar{h}),\dots ,\lambda_k(\bar{h})$ ，则必有唯一根，比如 $\lambda_{1}(\bar{h})\to 1$ (1是 $\rho (\lambda) = 0$ 的单根)当 $\bar{h}\rightarrow 0$ ，其余根和1保持一正距离；比如 $\left|\lambda_j(\bar{h}) - 1\right|\geqslant \delta_0 > 0,$ 当 $\bar{h}\to 0,j = 2,3,\dots ,k.$ 又 $u = \mathrm{e}^{\mu t}$ 满足 $u^{\prime} = \mu u$ .若方法为 $p$ 阶，则

$$
\sum_ {j = 0} ^ {k} \alpha_ {j} \mathrm {e} ^ {\mu (t + j h)} - \bar {h} \sum_ {j = 0} ^ {k} \beta_ {j} \mathrm {e} ^ {\mu (t + j h)} = O \left(h ^ {p + 1}\right) = O \left(\bar {h} ^ {p + 1}\right),
$$

从而

$$
\rho \left(\mathrm {e} ^ {\bar {h}}\right) - \bar {h} \sigma \left(\mathrm {e} ^ {\bar {h}}\right) = O \left(\bar {h} ^ {p + 1}\right).
$$

将左端分解因子，得

$$
\left(\mathrm {e} ^ {\bar {h}} - \lambda_ {1} (\bar {h})\right) \left(\mathrm {e} ^ {\bar {h}} - \lambda_ {2} (\bar {h})\right) \dots \left(\mathrm {e} ^ {\bar {h}} - \lambda_ {k} (\bar {h})\right) = O \left(\bar {h} ^ {p + 1}\right).
$$

因 $\exp (\bar{h})\to 1$ 当 $\bar{h}\rightarrow 0$ ，故除第一个因子外，其余因子均有正下界，即

$$
\left| \mathrm {e} ^ {\bar {h}} - \lambda_ {j} (\bar {h}) \right| \geqslant \delta_ {1} > 0, \quad j = 2, 3, \dots , k
$$

$(\delta_{1}$ 是与 $\bar{h}$ 无关的正数), 因此

$$
\lambda_ {1} (\bar {h}) = \mathrm {e} ^ {\bar {h}} + O \left(\bar {h} ^ {p + 1}\right) = 1 + \bar {h} + O \left(\bar {h} ^ {2}\right) + O \left(\bar {h} ^ {p + 1}\right). \tag {1.133}
$$

由方法的绝对稳定性知 $|\lambda_1(\bar{h})| < 1$ 恒成立, 又 $p \geqslant 1$ , 故必有 $\operatorname{Re}\bar{h} = \operatorname{Re}\mu h < 0$ , 因此 $\operatorname{Re}\mu < 0$ .

注意单位圆周 $|\lambda| = 1$ 在变换 (1.132) 之下, 所得的映像就是绝对稳定域 $D_A$ 的边界. 由 (1.133), $\bar{h} = 0$ 时 $\lambda_1(0) = 1$ , 故 $D_A$ 的边界经过 0.

检验绝对稳定性归结为检验特征方程（1.131）的根是否在单位圆内 $(|\lambda| < 1)$ ，这有很多判别法，如著名的 Routh-Hurwitz 准则，Schur 准则和 Miller 准则。这里只列一个简单而常用的判别法：复系数二次方程 $\lambda^2 - b\lambda - c = 0$ 的根在单位圆内的充要条件是

$$
2 - | b | ^ {2} - \left| b ^ {2} + 4 c \right| + 2 | c | ^ {2} > 0, \quad | c | <   1. \tag {1.134}
$$

# 1.5.3 应用例子

例1.5 Adams法的绝对稳定域

Adams外插法和内插法的特征多项式分别为

$$
\rho (\lambda) - \bar {h} \sigma (\lambda) = \lambda^ {k - 1} (\lambda - 1) - \bar {h} \sum_ {l = 0} ^ {k - 1} b _ {k - 1, l} \lambda^ {k - l - 1}
$$

和

$$
\rho (\lambda) - \bar {h} \sigma (\lambda) = \lambda^ {k - 1} (\lambda - 1) - \bar {h} \sum_ {l = 0} ^ {k} b _ {k l} ^ {*} \lambda^ {k - l},
$$

其中 $b_{kl}$ 和 $b_{kl}^{*}$ 分别由表1.2及表1.4给出. 可以算出 $k$ 步 $(k = 1,2,3,4)$ Adams外插法和内插法的绝对稳定域与实轴的交集分别为 $(\alpha_E,0)$ 和 $(\alpha_I,0)$ , 其中 $\alpha_E,\alpha_I$ 如表1.7.从表中看出，内插法的稳定域比外插法大

表 1.7 Adams 方法的绝对稳定域  

<table><tr><td>k</td><td>1</td><td>2</td><td>3</td><td>4</td></tr><tr><td>αE</td><td>-2</td><td>-1</td><td>-6/11</td><td>-3/10</td></tr><tr><td>αI</td><td>-∞</td><td>-6</td><td>-3</td><td>-90/49</td></tr></table>

例1.6 Runge-Kutta法的绝对稳定域

在 Runge-Kutta 法中取 $f = \mu u$ ( $\mu$ 为常数), 则

$$
k _ {1} = \mu u _ {n},
$$

$$
k _ {2} = \mu (1 + b _ {2 1} \mu h) u _ {n} = \mu P _ {1} (\mu h) u _ {n},
$$

$$
\begin{array}{l} k _ {3} = \mu \left(u _ {n} + h \sum_ {j = 1} ^ {2} b _ {3 j} k _ {j}\right) = \mu \left(1 + b _ {3 1} \mu h P _ {1} (\mu h)\right) u _ {n} = \mu P _ {2} (\mu h) u _ {n}, \\ k _ {m} = \mu P _ {m - 1} (\mu h) u _ {n}, \\ \end{array}
$$

其中 $P_{i}(\lambda)$ 是 $i$ 次多项式，从而

$$
\varphi \left(t _ {n}, u _ {n}, h\right) = \mu \left(\sum_ {i = 1} ^ {m} c _ {i} P _ {i - 1} (\mu h)\right) u _ {n}.
$$

Runge-Kutta法为

$$
u _ {n + 1} = u _ {n} + h \varphi (t _ {n}, u _ {n}, h) = u _ {n} + P _ {m} (\mu h) u _ {n}, \quad n = 0, 1, \dots . \tag {1.135}
$$

其次, 将 $u(t_{n+1})$ 在 $t_n$ 展开:

$$
u \left(t _ {n + 1}\right) = u \left(t _ {n}\right) + h u ^ {\prime} \left(t _ {n}\right) + \dots + \frac {h ^ {p}}{p !} u ^ {(p)} \left(t _ {n}\right) + O \left(h ^ {p + 1}\right).
$$

将 $u^{\prime} = \mu u$ 的解 $u(t) = \mathrm{e}^{\mu t}$ 代到上式得

$$
u \left(t _ {n + 1}\right) = \left(1 + \mu h + \dots + \frac {(\mu h) ^ {p}}{p !}\right) u \left(t _ {n}\right) + O \left(h ^ {p + 1}\right).
$$

若方法是 $p$ 阶的，则

$$
1 + P _ {m} (\mu h) = 1 + \mu h + \dots + \frac {(\mu h) ^ {p}}{p !},
$$

而 $\mu h$ 任意, 故 $m \geqslant p$ . 取 $m = p$ , 将 (1.135) 写为

$$
u _ {n + 1} = \left(\sum_ {l = 0} ^ {m} \frac {(\mu h) ^ {l}}{l !}\right) u _ {n}, \quad n = 0, 1, \dots ,
$$

其特征方程是一次的，唯一的特征值

$$
\lambda_ {1} = 1 + \bar {h} + \frac {1}{2 !} \bar {h} ^ {2} + \dots + \frac {1}{m !} \bar {h} ^ {m}, \quad \bar {h} = \mu h. \tag {1.136}
$$

注意, 当 $m = 1,2,3,4$ 时, 解不等式 $|\lambda_1| < 1$ 就可得出绝对稳定域 $D_A$ , 参看图1.2. 表1.8是 $D_A$ 与实轴相交的区间, 从中看出, Runge-Kutta法的绝对稳定域一般比线性多步法大, 这是它的优点.

在使用 Runge-Kutta 法时, 应取 $h$ 使 $\bar{h}$ 属于绝对稳定域, 否则有可能产生很大误差 (虽然方法稳定). 例如用四阶 Runge-Kutta 法 (1.119) 求解

$$
u ^ {\prime} = - 2 0 u, \quad u (0) = 1,
$$

步长为0.1及0.2.当步长取0.1时， $\bar{h}$ 属于绝对稳定域，误差随 $n$ 下降到零.当步长取0.2时， $\bar{h}$ 不属于绝对稳定域，误差很快增长.计算结果的误差如表1.9.

![](images/408cfd926f22009055c91803e12007192ca8029ac45506c7e4386bbda5a98510.jpg)  
图1.2 特征值的分布

表 1.8 Runge-Kutta 法的绝对稳定域  

<table><tr><td>级</td><td>λ1</td><td>绝对稳定域</td></tr><tr><td>一级</td><td>1+h</td><td>(-2,0)</td></tr><tr><td>二级</td><td>1+h+h2/2</td><td>(-2,0)</td></tr><tr><td>三级</td><td>1+h+h2/2+h3/6</td><td>(-2.51,0)</td></tr><tr><td>四级</td><td>1+h+h2/2+h3/6+h4/24</td><td>(-2.78,0)</td></tr></table>

表 1.9 四阶 Runge-Kutta 法的计算误差  

<table><tr><td>t</td><td>h=0.1</td><td>h=0.2</td></tr><tr><td>0.0</td><td>0</td><td>0</td></tr><tr><td>0.2</td><td>-0.092795</td><td>4.98</td></tr><tr><td>0.4</td><td>-0.012010</td><td>25.0</td></tr><tr><td>0.6</td><td>-0.001366</td><td>125.0</td></tr><tr><td>0.8</td><td>-0.000152</td><td>625.0</td></tr><tr><td>1.0</td><td>-0.000017</td><td>3125.0</td></tr></table>

Runge-Kutta法和预-校算法是求解常微分方程初值问题的有效数值方法.本节介绍的Runge-Kutta法都是显式的，为了进一步改善稳定域，人们也采用隐式Runge-Kutta法，参见[25,26,27].

# 1.5.4 习题

1. 证明方法 $u_{n + 1} - u_n = hf_{n + 1}$ 对所有 $\bar{h}\in (-\infty ,0)$ 绝对稳定  
2. 验证表 1.7 所列 Adams 外插法和内插法的绝对稳定域  
3. 求二级二阶隐式 Runge-Kutta 法

$$
\begin{array}{l} u _ {n + 1} = u _ {n} + \frac {1}{2} h (k _ {1} + k _ {2}), \\ k _ {1} = f \left(t _ {n}, u _ {n}\right), \\ k _ {2} = f \left(t _ {n} + h, u _ {n} + \frac {1}{2} h \left(k _ {1} + k _ {2}\right)\right) \\ \end{array}
$$

的绝对稳定域

# *1.6 一阶方程组和刚性问题

# 1.6.1 对一阶方程组的推广

实际中遇到的不只是一阶方程式, 还会有含 $m$ 个方程的一阶方程组, 其一般形式为

$$
\left\{ \begin{array}{l} \frac {\mathrm {d} u _ {1}}{\mathrm {d} t} = f _ {1} (t, u _ {1}, \dots , u _ {m}), \\ \frac {\mathrm {d} u _ {2}}{\mathrm {d} t} = f _ {2} (t, u _ {1}, \dots , u _ {m}), \\ \dots \dots \dots \dots \\ \frac {\mathrm {d} u _ {m}}{\mathrm {d} t} = f _ {m} (t, u _ {1}, \dots , u _ {m}), \end{array} \right. \tag {1.137}
$$

这里 $f_{1}, f_{2}, \dots, f_{m}$ ( $t_{0} \leqslant t \leqslant T, |u_{i}| < \infty, i = 1, 2, \dots, m$ ) 是 $m + 1$ 个变元的连续函数. 为使 (1.137) 的解确定, 还需给出初值条件:

$$
u _ {i} \left(t _ {0}\right) = u _ {0 i}, \quad i = 1, 2, \dots , m. \tag {1.138}
$$

这就是一阶方程组的初值问题. 如果实际问题不是一阶方程组而是高阶方程式, 我们也可

以把它化为一阶方程组，例如对于 $m$ 阶微分方程

$$
u ^ {(m)} = f (t, u, u ^ {\prime}, \dots , u ^ {(m - 1)}),
$$

只要引进新变量

$$
u _ {1} = u, \quad u _ {2} = u ^ {\prime}, \dots , \quad u _ {m} = u ^ {(m - 1)},
$$

就化为一阶方程组

$$
\left\{ \begin{array}{l} \frac {\mathrm {d} u _ {1}}{\mathrm {d} t} = u _ {2}, \\ \frac {\mathrm {d} u _ {2}}{\mathrm {d} t} = u _ {3}, \\ \dots \dots \dots \dots . \\ \frac {\mathrm {d} u _ {m}}{\mathrm {d} t} = f (t, u _ {1}, u _ {2}, \dots , u _ {m - 1}). \end{array} \right.
$$

此种转化不仅是理论上需要，在计算上也可能更方便.引进向量记号：

$$
\boldsymbol {u} = \left(u _ {1}, u _ {2}, \dots , u _ {m}\right) ^ {\mathrm {T}}, \quad \boldsymbol {f} = \left(f _ {1}, f _ {2}, \dots , f _ {m}\right) ^ {\mathrm {T}}, \quad \boldsymbol {u} _ {0} = \left(u _ {0 1}, u _ {0 2}, \dots , u _ {0 m}\right) ^ {\mathrm {T}}
$$

则(1.137)可写为向量形式：

$$
\boldsymbol {u} ^ {\prime} = \boldsymbol {f} (t, \boldsymbol {u}), \quad \boldsymbol {u} (t _ {0}) = \boldsymbol {u} _ {0}. \tag {1.139}
$$

若 $f(t,u)$ 关于 $\pmb{u}$ 满足Lipschitz条件,则问题(1.139)有唯一解(参见[2]).

前面介绍的线性多步法，预-校算法和Runge-Kutta法都可直接推广到一阶方程组，只需用向量代替相应的标量.所有关于相容性、稳定性和收敛性的定义和结论都可推广到方程组，只需将绝对值 $|\cdot |$ 换成 $m$ 维欧氏空间的向量模 $\| \cdot \|$ .例如解方程组的线性多步法是

$$
L \left[ \boldsymbol {u} _ {n}; h \right] = \sum_ {j = 0} ^ {k} \left[ \alpha_ {j} \boldsymbol {u} _ {n + j} - h \beta_ {j} \boldsymbol {f} \left(t _ {n + j}, \boldsymbol {u} _ {n + j}\right) \right] = 0, \tag {1.140}
$$

$\alpha_{j}, \beta_{j}$ 是标量, $L[u_{n}; h]$ 是向量算子. 稳定性的充要条件仍然是第一特征多项式满足根条件.

但绝对稳定性的定义要适当修改. 因为在方程式的情形, $\pmb{f}$ 的导数 $\frac{\partial f}{\partial u}$ 是标量, 在模型问题 $\pmb{u}^{\prime} = \mu \pmb{u}$ 中的 $\mu$ 也是标量. 现在 $\frac{\partial f}{\partial u}$ 是 $m$ 阶Jacobi矩阵, 应当用 $m$ 阶常矩阵 $A$ 代替 $\frac{\partial f}{\partial u}$ . 和 (1.126) 平行的线性化误差方程是

$$
\sum_ {j = 0} ^ {k} \left(\alpha_ {j} \mathbf {I} - h \beta_ {j} \mathbf {A}\right) \bar {\mathbf {e}} _ {n + j} = \boldsymbol {\eta} _ {n}. \tag {1.141}
$$

假定矩阵 $A$ 的特征值互异，则可对角化

$$
\boldsymbol {H} ^ {- 1} \boldsymbol {A} \boldsymbol {H} = \left[ \begin{array}{c c c c} \mu_ {1} & & & \\ & \mu_ {2} & & \\ & & \ddots & \\ & & & \mu_ {m} \end{array} \right].
$$

从而方程组 (1.141) 可化成 $m$ 个独立的方程式, 相应的特征多项式有 $m$ 个:

$$
\rho (\lambda) - h \mu_ {l} \sigma (\lambda) = 0, \quad l = 1, 2, \dots , m, \tag {1.142}
$$

其中

$$
\rho (\lambda) = \sum_ {j = 0} ^ {k} \alpha_ {j} \lambda^ {j}, \sigma (\lambda) = \sum_ {j = 0} ^ {k} \beta_ {j} \lambda^ {j}.
$$

以 $r_j^{(l)}(j = 1,2,\dots ,k,l = 1,2,\dots ,m)$ 表示（1.142）的根的绝对值，则可将方法的绝对稳定性定义为

$$
r _ {j} ^ {(l)} <   1, \quad j = 1, 2, \dots , k, \quad l = 1, 2, \dots , m. \tag {1.143}
$$

一般说来，矩阵 $\mathbf{A}$ 非对称，特征值 $\mu_{l}$ 是复的，所以由(1.143）确定的参数 $\bar{h} = h\mu$ 的稳定域也是复的.

# 1.6.2 刚性问题

一阶方程组的解在 $t\to \infty$ 的过程中，各分量的变化可能相差很大，这会给计算带来一定困难.例如考虑如下线性常微分方程组：

$$
\boldsymbol {u} ^ {\prime} = \boldsymbol {A} \boldsymbol {u}, \quad \boldsymbol {u} \in \mathbb {R} ^ {3}, \tag {1.144}
$$

$$
\boldsymbol {u} (0) = (2, 1, 2) ^ {\mathrm {T}}, \tag {1.145}
$$

其中

$$
\boldsymbol {A} = \left[ \begin{array}{c c c} - 0. 1 & - 4 9. 9 & 0 \\ 0 & - 5 0 & 0 \\ 0 & 7 0 & - 3 \times 1 0 ^ {4} \end{array} \right]. \tag {1.146}
$$

$\mathbf{A}$ 的特征值为

$$
\lambda_ {1} = - 0. 1, \quad \lambda_ {2} = - 5 0, \quad \lambda_ {3} = - 3 \times 1 0 ^ {4}.
$$

所以方程组（1.144）（1.145）的解为

$$
u _ {1} (t) = \mathrm {e} ^ {- 0. 1 t} + \mathrm {e} ^ {- 5 0 t},
$$

$$
u _ {2} (t) = \mathrm {e} ^ {- 5 0 t},
$$

$$
\begin{array}{l} u _ {3} (t) = \frac {7 0}{3 \times 1 0 ^ {4} - 7 0} \mathrm {e} ^ {- 5 0 t} + \left(2 - \frac {7 0}{3 \times 1 0 ^ {4} - 7 0}\right) \mathrm {e} ^ {- 3 \times 1 0 ^ {- 4}} \\ = 0. 0 0 2 3 3 7 \mathrm {e} ^ {- 5 0 t} + 1. 9 9 7 6 6 3 \mathrm {e} ^ {- 3 \times 1 0 ^ {4} t}. \\ \end{array}
$$

显然当 $t \to \infty$ 时解的各分量 $u_{i}(t)$ 均按指数衰减到稳态解 $u_{i}(t) = 0, i = 1,2,3$ . 但这些分量收敛到稳态解的速度很不一样, 这跟特征值的大小悬殊有关. 若用数值方法求解, 例如用 Runge-Kutta 法, 其绝对稳定域为 $(- \alpha_{R}, 0) = (-2.78, 0)$ . 为使方法绝对稳定, 应要求 $|h \lambda_{3}| = 3 \times 10^{4} h \leqslant \alpha_{R} = 2.78$ 或

$$
h \approx \frac {\alpha_ {R}}{| \lambda_ {3} |} = \frac {2 . 7 8}{3 \times 1 0 ^ {4}} \approx 1 0 ^ {- 4}.
$$

其次，为使 $u(t)$ 充分接近稳态解，应要求 $\mathrm{e}^{\lambda_1t} = \mathrm{e}^{-0.1t}$ 充分接近零，比如 $\mathrm{e}^{-0.1t}\approx \mathrm{e}^{-4}$ 或 $T\approx \frac{4}{|\lambda_1|}\approx 40,$ 于是计算步数

$$
N \approx \frac {T}{h} \approx \frac {4 | \lambda_ {3} |}{\alpha_ {R} | \lambda_ {1} |} \approx \frac {4 0}{1 0 ^ {- 4}} = 4 \times 1 0 ^ {5}.
$$

若注意到每步要计算右端函数 $f$ ，便知计算量是很大的

产生上述困难主要来自方程的系数矩阵（或一般方程的Jacobi矩阵）的特征值的分布，它们位于左半平面，而按绝对值相差悬殊，这类方程称为刚性(stiff)方程.

考虑线性常系数系统

$$
\boldsymbol {u} ^ {\prime} (t) = \boldsymbol {A} \boldsymbol {u} + \boldsymbol {g} (t) \tag {1.147}
$$

和非线性系统

$$
\boldsymbol {u} ^ {\prime} (t) = \boldsymbol {f} (t, \boldsymbol {u}) \tag {1.148}
$$

定义1.5 称(1.147)为刚性的，如果 $\lambda_{i}$ 是矩阵 $\mathbf{A}$ 的特征值，满足

(i) $\operatorname{Re}\lambda_i < 0, i = 1,2,\dots,m,$   
$\max_{i}|\operatorname {Re}\lambda_{i}|$ (ii) $\frac{\min_{i}|\operatorname{Re}\lambda_i|}{\min_{i}} = R\gg 1,$ 其中 $\operatorname {Re}\lambda_{i}$ 是 $\lambda_{i}$ 的实部， $R$ 称为刚性比

定义1.6 称(1.148)为刚性的，如果在 $t$ 的区间 $I = [0,T]$ 上， $\pmb{f}$ 的Jacobi矩阵 $\frac{\partial f}{\partial u}$ 的特征值 $\lambda_{i}(t)$ 满足定义1.5的条件(i)(ii).

粗略说来，刚性方程组的特征是其解同时存在快变和慢变部分，这类方程组在生物学、化学、电子学和控制论等领域有重要应用.

# 1.6.3 A稳定性

则为保证方法绝对稳定，最好要求绝对稳定域就是左复平面 $Z_{-}$ ： $\operatorname {Re}(h\lambda) <   0$ ，为此引进

A 稳定定义.

定义1.7 线性多步法说是A稳定的，如果将它用于模型问题

$$
u ^ {\prime} (t) = \lambda u (t) \tag {1.149}
$$

的绝对稳定域就是左复平面 $Z_{-}$ ，其中 $\lambda$ 是复数

显然用A稳定算法计算时，步长 $h$ 不必受绝对稳定条件的限制

为了判别线性多步法A稳定，将它用于模型问题(1.149)，得线性差分方程：

$$
\sum_ {j = 0} ^ {k} \left(\alpha_ {j} - \bar {h} \beta_ {j}\right) u _ {n + j} = 0, \quad \bar {h} = h \lambda \tag {1.150}
$$

相应的特征方程为 $\rho (\lambda) - \bar{h}\sigma (\lambda) = 0$ 或

$$
\bar {h} = \frac {\rho (\lambda)}{\sigma (\lambda)}, \tag {1.151}
$$

其中

$$
\rho (\lambda) = \sum_ {j = 0} ^ {k} \alpha_ {j} \lambda^ {j}, \quad \sigma (\lambda) = \sum_ {j = 0} ^ {m} \beta_ {j} \lambda^ {j}, \quad m \leqslant k. \tag {1.152}
$$

(1.151) 是一个由 $\lambda \in Z$ 到 $\bar{h} \in Z$ 的单值解析变换, 其逆是由 $\bar{h} \in Z$ 到 $\lambda \in Z$ 的 $k$ 值逆变换. 由定义推出

命题1.2 设 $\lambda_{i}(i = 1,2,\dots ,k)$ 是方程(1.151）的根，则下列表述等价：

(i) 线性多步法A稳定；  
(ii) $\operatorname{Re}\bar{h} < 0 \Rightarrow |\lambda_i| < 1, i = 1,2,\dots,k;$   
(iii) $|\lambda |\geqslant 1\Rightarrow \operatorname {Re}\bar{h} (\lambda)\geqslant 0.$

例1.7 Euler向后公式

$$
u _ {n + 1} = u _ {n} + h f _ {n + 1}.
$$

由于 $\rho (\lambda) = \lambda -1,\sigma (\lambda) = \lambda$ ，则

$$
\operatorname {R e} \bar {h} (\lambda) = \operatorname {R e} \frac {\lambda - 1}{\lambda} = \frac {| \lambda | ^ {2} - | \lambda | \cos \theta}{| \lambda | ^ {2}} = \frac {| \lambda | (| \lambda | - \cos \theta)}{| \lambda | ^ {2}}.
$$

显然当 $|\lambda |\geqslant 1$ 时， $\operatorname {Re}\bar{h} (\lambda)\geqslant 0,$ 故Euler向后公式A稳定

例1.8 梯形公式

$$
u _ {n + 1} = u _ {n} + \frac {h}{2} \left(f _ {n + 1} + f _ {n}\right).
$$

因 $\rho (\lambda) = \lambda -1,\sigma (\lambda) = \frac{1}{2} (\lambda +1)$ ，故

$$
\operatorname {R e} \bar {h} (\lambda) = \operatorname {R e} \frac {\rho (\lambda)}{\sigma (\lambda)} = \operatorname {R e} \frac {\lambda - 1}{(\lambda + 1) / 2} = 2 \left(\frac {| \lambda | ^ {2} - 1}{| \lambda + 1 | ^ {2}}\right).
$$

于是当 $|\lambda |\geqslant 1$ 时， $\operatorname {Re}\bar{h} (\lambda)\geqslant 0.$ 所以梯形公式A稳定

例1.9 考虑 $k$ 步线性方法：

$$
u _ {n + k} - u _ {n} = \frac {h}{2} k \left(f _ {n + k} + f _ {n}\right).
$$

因

$$
\rho (\lambda) = \lambda^ {k} - 1, \quad \sigma (\lambda) = \frac {k}{2} \left(\lambda^ {k} + 1\right),
$$

故

$$
\bar {h} = \frac {\rho (\lambda)}{\sigma (\lambda)} = \frac {2}{k} \frac {\lambda^ {k} - 1}{\lambda^ {k} + 1} = \frac {2}{k} \frac {| \lambda | ^ {2 k} - 1 + \mathrm {i} 2 | \lambda | ^ {k} \sin \theta}{| \lambda^ {k} + 1 | ^ {2}}.
$$

于是

$$
\operatorname {R e} \bar {h} (\lambda) = \frac {2}{k} \frac {| \lambda | ^ {2 k} - 1}{| \lambda^ {k} + 1 | ^ {2}}.
$$

显然当 $|\lambda |\geqslant 1$ 时， $\operatorname {Re}\bar{h} (\lambda)\geqslant 0,$ 故此 $k$ 步法A稳定

现在看Euler法

$$
u _ {n + 1} - u _ {n} = h f _ {n}.
$$

这是显方法, 相应的 $\rho (\lambda) = \lambda -1,\sigma (\lambda) = 1$ ，此时

$$
\operatorname {R e} \bar {h} = \operatorname {R e} \frac {1}{\lambda - 1} = \operatorname {R e} \frac {\operatorname {R e} \lambda - 1 + \mathrm {i} \operatorname {I m} \lambda}{| \lambda - 1 | ^ {2}} = \frac {\operatorname {R e} \lambda - 1}{| \lambda - 1 | ^ {2}}.
$$

显然若 $|\lambda| \geqslant 1$ , 而 $\operatorname{Re} \lambda < 1$ , 则 $\operatorname{Re} \bar{h} < 0$ , 所以 Euler 法非 A 稳定.

实际上，可以证明显线性多步法都不是A稳定的(参见[9]).

# 1.6.4 数值例子

1. 用向后 Euler 公式求解初值问题 (1.144) — (1.146). 以 $u_{j}^{n}$ 表示 $u_{j}(t)$ 于 $t_{n} = nh$ 的近似, 则计算公式为

$$
u _ {1} ^ {n + 1} - u _ {1} ^ {n} = h \left(- 0. 1 u _ {1} ^ {n + 1} - 4 9. 9 u _ {2} ^ {n + 1}\right),
$$

$$
u _ {2} ^ {n + 1} - u _ {2} ^ {n} = h (- 5 0 u _ {2} ^ {n + 1}),
$$

$$
u _ {3} ^ {n + 1} - u _ {3} ^ {n} = h \left(7 0 u _ {2} ^ {n + 1} - 3 \times 1 0 ^ {4} u _ {3} ^ {n + 1}\right)
$$

或

$$
(1 + 0. 1 h) u _ {1} ^ {n + 1} + 4. 9 9 h u _ {2} ^ {n + 1} = u _ {1} ^ {n},
$$

$$
(1 + 5 0 h) u _ {2} ^ {n + 1} = u _ {2} ^ {n},
$$

$$
\left(1 + 3 \times 1 0 ^ {4} h\right) u _ {3} ^ {n + 1} - 7 0 h u _ {2} ^ {n + 1} = u _ {3} ^ {n}.
$$

以 $h = 1$ 和 $n = 1,2,\dots$ 计算，直至

$$
\| u ^ {n} \| = \left[ \left(u _ {1} ^ {n}\right) ^ {2} + \left(u _ {2} ^ {n}\right) ^ {2} + \left(u _ {3} ^ {n}\right) ^ {2} \right] ^ {\frac {1}{2}} \leqslant \mathrm {e} ^ {- 4} \approx 0. 0 0 0 0 2 7.
$$

计算结果如表1.10.从表中看到，当 $t = 100$ 时达到稳定解

表 1.10 向后 Euler 公式的计算结果  

<table><tr><td rowspan="2">t</td><td colspan="3">精确解</td><td colspan="3">数值解</td></tr><tr><td>u1(t)</td><td>u2(t)</td><td>u3(t)</td><td>u1(t)</td><td>u2(t)</td><td>u3(t)</td></tr><tr><td>0</td><td>2.000 000</td><td>1.000 000</td><td>2.000 000</td><td>2.000 000</td><td>1.000 000</td><td>2.000 000</td></tr><tr><td>1</td><td>0.904 837</td><td>0.000 000</td><td>0.000 000</td><td>0.928 698</td><td>0.019 607</td><td>0.000 112</td></tr><tr><td>2</td><td>0.818 730</td><td>0.000 000</td><td>0.000 000</td><td>0.826 830</td><td>0.000 384</td><td>0.000 000</td></tr><tr><td>3</td><td>0.740 818</td><td>0.000 000</td><td>0.000 000</td><td>0.751 322</td><td>0.000 007</td><td>0.000 000</td></tr><tr><td>4</td><td>0.670 320</td><td>0.000 000</td><td>0.000 000</td><td>0.683 013</td><td>0.000 000</td><td>0.000 000</td></tr><tr><td>5</td><td>0.606 530</td><td>0.000 000</td><td>0.000 000</td><td>0.620 921</td><td>0.000 000</td><td>0.000 000</td></tr><tr><td>:</td><td>:</td><td>:</td><td>:</td><td>:</td><td>:</td><td>:</td></tr><tr><td>95</td><td>0.000 082</td><td>0.000 000</td><td>0.000 000</td><td>0.000 128</td><td>0.000 000</td><td>0.000 000</td></tr><tr><td>96</td><td>0.000 074</td><td>0.000 000</td><td>0.000 000</td><td>0.000 116</td><td>0.000 000</td><td>0.000 000</td></tr><tr><td>97</td><td>0.000 067</td><td>0.000 000</td><td>0.000 000</td><td>0.000 106</td><td>0.000 000</td><td>0.000 000</td></tr><tr><td>98</td><td>0.000 061</td><td>0.000 000</td><td>0.000 000</td><td>0.000 096</td><td>0.000 000</td><td>0.000 000</td></tr><tr><td>99</td><td>0.000 055</td><td>0.000 000</td><td>0.000 000</td><td>0.000 087</td><td>0.000 000</td><td>0.000 000</td></tr><tr><td>100</td><td>0.000 050</td><td>0.000 000</td><td>0.000 000</td><td>0.000 087</td><td>0.000 000</td><td>0.000 000</td></tr></table>

2. 用向前 Euler 公式求解初值问题 (1.144) — (1.146), 计算公式为

$$
u _ {1} ^ {n + 1} - u _ {1} ^ {n} = - h \left(0. 1 h u _ {1} ^ {n} + 4. 9 9 u _ {2} ^ {n}\right),
$$

$$
u _ {2} ^ {n + 1} - u _ {2} ^ {n} = - 5 0 h u _ {2} ^ {n},
$$

$$
u _ {3} ^ {n + 1} - u _ {3} ^ {n} = h \left(7 0 u _ {2} ^ {n} - 3 \times 1 0 ^ {4} u _ {3} ^ {n}\right).
$$

取 $h = 0.5$ ，从表1.11看出，计算到 $t = 2$ 即出现不稳定现象

表 1.11 向前 Euler 公式的计算结果  

<table><tr><td rowspan="2">t</td><td colspan="3">精确解</td><td colspan="3">数值解</td></tr><tr><td>u1(t)</td><td>u2(t)</td><td>u3(t)</td><td>u1(t)</td><td>u2(t)</td><td>u3(t)</td></tr><tr><td>0</td><td>2.000 000</td><td>1.000 000</td><td>2.000 000</td><td>2.000 000</td><td>1.000 000</td><td>2.000 000</td></tr><tr><td>0.5</td><td>0.951 229</td><td>0.000 000</td><td>0.000 000</td><td>-23.050 000</td><td>-24.000 000</td><td>-29.963 000</td></tr><tr><td>1</td><td>0.904 837</td><td>0.000 000</td><td>0.000 000</td><td>576.902 499</td><td>576.000 000</td><td>4.494 141 × 10^8</td></tr><tr><td>1.5</td><td>0.860 707</td><td>0.000 000</td><td>0.000 000</td><td>-13 823.142 625</td><td>-13 824.000 000</td><td>-6.740 763 × 10^12</td></tr><tr><td>2</td><td>0.818 730</td><td>0.000 000</td><td>0.000 000</td><td>3.317 768 × 10^5</td><td>3.317 760 × 10^5</td><td>1.011 047 × 10^17</td></tr></table>

# *1.7 外推法

![](images/fa1c7da3a0c53b1969defc55b5bc81cd9d83f89d11ea7bb026419ac898b4ce2a.jpg)

人物简介

# 1.7.1 多项式外推

设 $h > 0$ 是离散化参数 (比如是步长), $A(h)$ 是按某种算法得到的数 $A_0$ 的近似, 且 $\lim_{h\to 0}A(h) = A_0$ (定义 $A(0) = A_0$ ). 假定对任意正整数 $N$ , $A(h)$ 有渐近展开式

$$
A (h) = A _ {0} + A _ {1} h + A _ {2} h ^ {2} + \dots + A _ {N} h ^ {N} + O \left(h ^ {N + 1}\right), \quad h \rightarrow 0, \tag {1.153}
$$

其中系数 $A_0, A_1, \dots, A_N$ 与 $h$ 无关. 若已算出 $A(h_0)$ 和 $A\left(\frac{1}{2} h_0\right)$ , 则由 (1.7.1) 知

$$
2 A \left(\frac {1}{2} h _ {0}\right) - A \left(h _ {0}\right) = A _ {0} + O \left(h _ {0} ^ {2}\right). \tag {1.154}
$$

可见左端较 $A(h_0)$ 和 $A\left(\frac{1}{2} h_0\right)$ 逼近 $A_0$ 的阶更高. 这种由已知近似组合成更好近似的方法称为外推法 (extrapolation methods), 其基本思想属于Richardson (1927). 计算数值积分的Romberg 方法实际上也是一种外推法. 如果除 $A(h_0)$ , $A\left(\frac{1}{2} h_0\right)$ 外还算出 $A\left(\frac{1}{4} h_0\right)$ , 则可找到有逼近阶更高的线性组合:

$$
\frac {1}{3} A \left(h _ {0}\right) - 2 A \left(\frac {h _ {0}}{2}\right) + \frac {8}{3} A \left(\frac {h _ {0}}{4}\right) = A _ {0} + O \left(h _ {0} ^ {3}\right). \tag {1.155}
$$

将这一思想推广，就得到更一般的所谓多项式外推法. 考虑一般的 $h$ 序列：

$$
h _ {0} > h _ {1} > h _ {2} > \dots > h _ {J} > 0.
$$

若展开式(1.153)成立，则总可求得 $A(h_j)(j = 0,1,\dots ,J)$ 的线性组合，使

$$
\sum_ {j = 0} ^ {J} c _ {j J} A \left(h _ {j}\right) = A _ {0} + O \left(h _ {0} ^ {J + 1}\right). \tag {1.156}
$$

实际上，求一以 $(h_j,A(h_j))(j = 0,1,\dots ,J)$ 为型值的 $A(h)$ 的 $J$ 次插值多项式 $P_{J}(h)$ 因

$$
A (h) = P _ {J} (h) + O \left(h _ {0} ^ {J + 1}\right),
$$

$A(0) = A_0, P_J(0)$ 是 $A(h_j)$ 的一次组合，于上式中令 $h = 0$ 即得 (1.156). 求 $P_J(h)$ 可用 Aitken 的逐步线性插值实现. 例如以 $\left(h_0, a_0^{(0)}\right)$ 和 $\left(h_1, a_1^{(0)}\right)$ 为型值的 $h$ 的一次多项式为

$$
I _ {0 1} (h) = \frac {1}{h _ {1} - h _ {0}} \left| \begin{array}{c c} a _ {0} ^ {(0)} & h _ {0} - h \\ a _ {1} ^ {(0)} & h _ {1} - h \end{array} \right|,
$$

其中 $a_0^{(0)} = A(h_0), a_1^{(0)} = A(h_1)$ . 令 $a_0^{(1)} = I_{01}(0)$ ，则由（1.156）知 $a_0^{(1)} = A_0 + O(h_0^2)$ . 特别地，当 $h_1 = \frac{1}{2} h_0$ 时， $a_0^{(1)}$ 就是（1.154）的左端. 同样过型值 $\left(h_1, a_1^{(0)}\right)$ ， $\left(h_2, a_2^{(0)}\right)$ $\left(a_2^{(0)} = A(h_2)\right)$ 的一次多项式为

$$
I _ {1 2} (h) = \frac {1}{h _ {2} - h _ {1}} \left| \begin{array}{c c} a _ {1} ^ {(0)} & h _ {1} - h \\ a _ {2} ^ {(0)} & h _ {2} - h \end{array} \right|.
$$

令 $a_1^{(1)} = I_{12}(0)$ 以 $I_{012}(h)$ 表示过三个点 $\left(h_0,a_0^{(0)}\right),\left(h_1,a_1^{(0)}\right),\left(h_2,a_2^{(0)}\right)$ 的二次多项式,则

$$
I _ {0 1 2} (h) = \frac {1}{h _ {2} - h _ {0}} \left| \begin{array}{c c} I _ {0 1} (h) & h _ {0} - h \\ I _ {1 2} (h) & h _ {2} - h \end{array} \right|.
$$

令 $a_0^{(2)} = I_{12}(0)$ ，则由（1.156）知 $a_0^{(2)} = A_0 + O\left(h_0^3\right)$ .特别地，当 $h_1 = \frac{h_0}{2},h_2 = \frac{h_0}{4}$ 时， $a_0^{(2)}$ 就是（1.155）的左端.如此可继续作下去，但插值多项式的次数不能超过方法的阶.实际应用时，作一两次插值就够了.

若渐近展开式形如

$$
A (h) = A _ {0} + A _ {2} h ^ {2} + A _ {4} h ^ {4} + O \left(h ^ {6}\right),
$$

则线性组合

$$
\frac {4}{3} A \left(\frac {1}{2} h _ {0}\right) - \frac {1}{3} A \left(h _ {0}\right) = A _ {0} + O \left(h _ {0} ^ {4}\right) \tag {1.157}
$$

比(1.154)的收敛阶更高, 可见外推法的效果和渐近展开的类型有关.

# 1.7.2 对初值问题的应用

设 $u(t)$ 是初值问题

$$
u ^ {\prime} = f (t, u), \quad u \left(t _ {0}\right) = u _ {0}
$$

的解, $u(t,h)$ 是由某一数值方法 (如线性多步法, 预-校算法, Runge-Kutta 法等) 确定的 $u(t)$ 的近似解, $h$ 是步长. 假定 $H$ 是基本步长 (可取大一些), $h_0 = \frac{H}{N_0}$ ( $N_0 \geqslant 1$ ). 从 $t = 0$ 出发, 用数值方法算 $N_0$ 步, 得到 $u(t_0 + H, h_0)$ . 然后取 $h_1 = \frac{H}{N_1}$ , $N_1 > N_0$ , 从 $t = t_0$ 出发用数值方法算 $N_1$ 步, 得到近似解 $u(t_0 + H, h_1)$ . 一般取 $h_j = \frac{H}{N_j}$ , $j = 0,1,\dots,J$ , $N_0 < N_1 < \dots < N_J$ , 从 $t = t_0$ 出发用数值方法算 $N_j$ 步, 得到近似解 $u(t_0 + H, h_j)$ . 如

果数值解有渐近展开：

$$
u (t, h) = u (t) + A _ {1} h + A _ {2} h ^ {2} + A _ {3} h ^ {3} + \dots + A _ {N} h ^ {N} + O \left(h ^ {N + 1}\right), \quad h \rightarrow 0, \tag {1.158}
$$

令 $a_{j}^{(0)} = u(t_{0} + H,h_{j}),j = 0,1,\dots ,J,$ 则可按逐步线性插值法得到 $u(t)$ 的更为精确的近似.通常取 $h_0 = H$ ， $h_1 = \frac{1}{2} h_0$ ，则

$$
\begin{array}{l} u \left(H, h _ {0}\right) = u (H) + O \left(h _ {0}\right), \\ 2 u (H, h _ {1}) - u (H, h _ {0}) = u (H) + O \left(h _ {0} ^ {2}\right). \tag {1.159} \\ \end{array}
$$

对单步法, 常有展开式

$$
u (t, h) = u (t) + \varepsilon (t) h ^ {p} + O \left(h ^ {p + 1}\right), \tag {1.160}
$$

其中 $\varepsilon (t)$ 为误差主项系数.用 $\frac{1}{2} h$ 为步长由 $t_0$ 算到 $t,$ 则

$$
u \left(t, \frac {1}{2} h\right) = u (t) + \varepsilon (t) \left(\frac {h}{2}\right) ^ {p} + O \left(h ^ {p + 1}\right). \tag {1.161}
$$

与 (1.160) 作线性组合 (线性插值), 得

$$
\frac {2 ^ {p} u \left(t , \frac {h}{2}\right) - u (t , h)}{2 ^ {p} - 1} = u (t) + O \left(h ^ {p + 1}\right). \tag {1.162}
$$

# 1.7.3 用外推法估计误差

假定方法的数值解 $u(t,h)$ 有渐进展开式(1.160).以 $e(t,h) = u(t) - u(t,h)$ 表示整体误差.由(1.160)，(1.161）得

$$
e \left(t, \frac {h}{2}\right) = 2 ^ {- p} e (t, h) + O \left(h ^ {p + 1}\right).
$$

而 $e\left(t,\frac{h}{2}\right) = u(t) - u\left(t,\frac{h}{2}\right) = e(t,h) + u(t,h) - u\left(t,\frac{h}{2}\right)$ , 代入上式左端, 解出

$$
e (t, h) = \frac {2 ^ {p}}{2 ^ {p} - 1} \left[ u \left(t, \frac {h}{2}\right) - u (t, h) \right] + O \left(h ^ {p + 1}\right). \tag {1.163}
$$

取

$$
\bar {e} (t, h) = \frac {2 ^ {p}}{2 ^ {p} - 1} \left[ u \left(t, \frac {h}{2}\right) - u (t, h) \right] \tag {1.164}
$$

作为 $e(t,h)$ 的近似，右端是可以用近似解估算的.利用(1.164)可以估计解的误差

例1.10 用Euler法解初值问题：

$$
u ^ {\prime} = - u, \quad u (0) = 1.
$$

其精确解 $u(t) = \mathrm{e}^{-t}$ . Euler 法的误差阶为 $O(h)$ , 外推公式 (1.159) 的误差阶为 $O(h^{2})$ . 表 1.12 就 $h = 2^{-j}$ ( $j = 0, 1, 2, 3, 4, 5, 6$ ) 列出 $t = 1$ 的近似值和误差. 从表 1.12 中看出, 外推解比 Euler 解精确, 用 (1.164) 作出的误差估计与实际误差基本符合.

表 1.12 $t = 1$ 的近似值和误差  

<table><tr><td>h</td><td>Euler 解 u(1,h)</td><td>2u(1,1/2)-u(1,h)</td><td>Euler 解的误差</td><td>外推解误差</td><td>用 (1.164) 估计出的误差</td></tr><tr><td>1</td><td>0.000000</td><td>0.500000</td><td>-0.367879</td><td>0.132121</td><td>-0.500000</td></tr><tr><td>1/2</td><td>0.250000</td><td>0.382813</td><td>-0.117879</td><td>0.014933</td><td>-0.132813</td></tr><tr><td>1/4</td><td>0.316406</td><td>0.370812</td><td>-0.051473</td><td>0.002932</td><td>-0.054405</td></tr><tr><td>1/8</td><td>0.343609</td><td>0.368539</td><td>-0.024271</td><td>0.000660</td><td>-0.024930</td></tr><tr><td>1/16</td><td>0.356074</td><td>0.368036</td><td>-0.011805</td><td>0.000157</td><td>-0.011962</td></tr><tr><td>1/32</td><td>0.362055</td><td>0.367918</td><td>-0.005824</td><td>0.000038</td><td>-0.005863</td></tr></table>

# 1.7.4 习题

1. 证明任一相容的显式单步法的解 $u(t, h)$ 有渐近展开式 (1.160).  
2. 设单位圆内正 $n$ 边型的半周长为 $A(h), nh = 1$ , 证明 $A(h)$ 可展成

$$
A (h) = \pi + A _ {2} h ^ {2} + A _ {4} h ^ {4} + \dots .
$$

并利用 $A(h), h = \frac{1}{4}, \frac{1}{6}, \frac{1}{8}$ 推出 $\pi$ 的近似式

<table><tr><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr><td></td><td></td><td></td><td></td><td></td><td></td></tr><tr><td></td><td></td><td></td><td></td><td></td><td></td></tr></table>

# 第2章

# 椭圆型方程的有限差分法

从本章开始，我们介绍偏微分方程的数值解法，主要是有限差分法和 Galerkin 有限元法。由于计算机只能存储有限个数据和作有限次运算，所以任何一种用计算机解题的方法，都必须把连续问题（微分方程的边值问题、初值问题等）离散化，最终化成有限形式的线性代数方程组。用差分法和有限元法将连续问题离散化的步骤是，首先对求解区域作网格剖分，用有限个网格节点代替连续区域，然后将微分算子离散化，从而把微分方程的定解问题化为线性代数方程组的求解问题。差分法和有限元法的主要区别是离散化的第二步。前者从定解问题的微分形式或积分形式出发，用数值微商或数值积分公式导出相应的线性代数方程组。后者从定解问题的变分形式出发，用 Ritz-Galerkin 法导出相应的线性代数方程组，但基函数按特定方式选取。本章和第 3 章、第 4 章讨论有限差分法，第 5 章、第 6 章讨论 Galerkin 有限元法。

差分法的基本问题有：

(1) 对求解域作网格剖分.   
一维情形是把区间分成一些等距或不等距的小区间，称之为单元。二维情形则把区域分割成一些均匀或不均匀的矩形，其边与坐标轴平行。也可分割成一些小三角形或凸四边形等。  
(2) 构造逼近微分方程定解问题的差分格式  
我们将介绍两种构造差分格式的方法：直接差分化法和有限体积法(积分插值法).  
(3) 差分解的存在唯一性、收敛性及稳定性的研究  
(4) 差分方程的解法

# 2.1 差分逼近的基本概念

考虑二阶常微分方程边值问题：

$$
L u = - \frac {\mathrm {d} ^ {2} u}{\mathrm {d} x ^ {2}} + q u = f, \quad a <   x <   b, \tag {2.1}
$$

$$
u (a) = \alpha , \quad u (b) = \beta , \tag {2.2}
$$

其中 $q, f$ 为 $[a, b]$ 上的连续函数, $q \geqslant 0$ ; $\alpha, \beta$ 为给定常数. 这是最简单的椭圆型方程第一边值问题.

将区间 $[a,b]$ 分成 $N$ 等分，分点为

$$
x _ {i} = a + i h, \quad i = 0, 1, \dots , N,
$$

$h = \frac{b - a}{N}$ . 于是我们得到区间 $I = [a, b]$ 的一个网格剖分. $x_{i}$ 称为网格的节点, $h$ 称为步长.

现在将方程 (2.1) 在节点 $x_{i}$ 离散化. 为此, 对充分光滑的解 $u$ , 由 Taylor 展开式可得

$$
\frac {u \left(x _ {i + 1}\right) - 2 u \left(x _ {i}\right) + u \left(x _ {i - 1}\right)}{h ^ {2}} = \left[ \frac {\mathrm {d} ^ {2} u (x)}{\mathrm {d} x ^ {2}} \right] _ {i} + \frac {h ^ {2}}{1 2} \left[ \frac {\mathrm {d} ^ {4} u (x)}{\mathrm {d} x ^ {4}} \right] _ {i} + O \left(h ^ {3}\right), \tag {2.3}
$$

其中 $[\cdot ]_i$ 表示方括号内的函数在 $x_{i}$ 点取值.于是在 $x_{i}$ 可将方程(2.1)写成

$$
- \frac {u \left(x _ {i + 1}\right) - 2 u \left(x _ {i}\right) + u \left(x _ {i - 1}\right)}{h ^ {2}} + q \left(x _ {i}\right) u \left(x _ {i}\right) = f \left(x _ {i}\right) + R _ {i} (u), \tag {2.4}
$$

其中

$$
R _ {i} (u) = - \frac {h ^ {2}}{1 2} \left[ \frac {\mathrm {d} ^ {4} u (x)}{\mathrm {d} x ^ {4}} \right] _ {i} + O \left(h ^ {3}\right). \tag {2.5}
$$

显然，当 $h\to 0$ 时， $R_{i}(u)$ 是 $h$ 的二阶无穷小量.若舍去 $R_{i}(u)$ ，则得到逼近方程(2.1）的差分方程：

$$
L _ {h} u _ {i} = - \frac {u _ {i + 1} - 2 u _ {i} + u _ {i - 1}}{h ^ {2}} + q _ {i} u _ {i} = f _ {i}, \tag {2.6}
$$

式中 $q_{i} = q(x_{i}),f_{i} = f(x_{i})$ .称 $R_{i}(u)$ 为差分方程(2.6）的截断误差.利用差分算子 $L_{h}$ 可将（2.4）写成

$$
L _ {h} u \left(x _ {i}\right) = f \left(x _ {i}\right) + R _ {i} (u). \tag {2.7}
$$

而在节点 $x_{i}$ 处，微分方程（2.1）为

$$
[ L u ] _ {i} = f (x _ {i}),
$$

与 (2.7) 相减, 得

$$
R _ {i} (u) = L _ {h} u \left(x _ {i}\right) - [ L u ] _ {i}. \tag {2.8}
$$

所以 $R_{i}(u)$ 是用差分算子 $L_{h}$ 逼近微分算子 $L$ 所引起的截断误差，在这里关于 $h$ 的阶为 $O(h^{2})$

差分方程 (2.6) 当 $i = 1,2,\dots ,N - 1$ 时成立，加上边值条件 $u_{0} = \alpha ,u_{N} = \beta$ ，就得到关于 $u_{i}$ 的线性代数方程组：

$$
\begin{array}{l} L _ {h} u _ {i} = - \frac {u _ {i + 1} - 2 u _ {i} + u _ {i - 1}}{h ^ {2}} + q _ {i} u _ {i} = f _ {i}, \quad i = 1, 2, \dots , N - 1, (2.9) \\ u _ {0} = \alpha , \quad u _ {N} = \beta . (2.10) \\ \end{array}
$$

它的解 $u_{i}$ 是 $u(x)$ 于 $x = x_{i}$ 的近似.称(2.9)—(2.10）为逼近(2.1)—(2.2）的差分方程或差分格式.由于（2.9）是用二阶中心差商代替（2.1）中二阶微商得到的，所以也称(2.9)—(2.10）为中心差分格式

注意方程 (2.9) 的个数等于网格内点 $x_{1}, x_{2}, \cdots, x_{N-1}$ 的个数 $N-1$ , 因此它是 $N-1$ 阶的方程组. 一般说来, 这是高阶方程组, 例如取 $N=100$ (即把区间 $[a, b]$ 作 100 等分), 则阶数为 99. 但每个方程的未知数最多有三个, 因此系数矩阵 $\mathbf{A}$ 的大量元素是零. 如

如果我们把方程和未知数按由左到右的节点顺序排列, 则 $\mathbf{A}$ 是对称的三对角矩阵. 例如取 $N = 5$ , 则

$$
\boldsymbol {A} = \left[ \begin{array}{c c c c} \frac {2}{h ^ {2}} + q _ {1} & - \frac {1}{h ^ {2}} & 0 & 0 \\ - \frac {1}{h ^ {2}} & \frac {2}{h ^ {2}} + q _ {2} & - \frac {1}{h ^ {2}} & 0 \\ 0 & - \frac {1}{h ^ {2}} & \frac {2}{h ^ {2}} + q _ {3} & - \frac {1}{h ^ {2}} \\ 0 & 0 & - \frac {1}{h ^ {2}} & \frac {2}{h ^ {2}} + q _ {4} \end{array} \right]
$$

我们可用消元法或迭代法求解方程组 (2.9) (2.10).

对于差分方程 (2.9) (2.10), 我们自然关心它是否有唯一解. 其次, 差分解 $u_{i}$ 当网格无限加密, 或者说 $h \to 0$ 时是否收敛到精确解 $u(x_{i})$ , 以及在何种度量意义下收敛, 收敛的速度如何. 为此需要引进若干记号.

以 $I_h$ 表示网格内点 $x_1, x_2, \dots, x_{N-1}$ 的集合， $\bar{I}_h$ 表示网格内点和界点 $x_0 = a, x_N = b$ 的集合。定义在 $I_h$ （相应的 $\bar{I}_h$ ）上的函数 $u_h(x_i) = u_i$ 称为 $I_h$ （相应的 $\bar{I}_h$ ）上的网格函数。和连续变量的函数类似，我们对 $I_h$ 上的网格函数引进范数

$$
\left\| u _ {h} \right\| _ {C} = \max  _ {1 \leqslant i \leqslant N - 1} | u _ {i} |, \tag {2.11}
$$

$$
\left\| u _ {h} \right\| _ {0} ^ {2} = \sum_ {i = 1} ^ {N - 1} h u _ {i} ^ {2}. \tag {2.12}
$$

$$
\left\| u _ {h} \right\| _ {1} ^ {2} = \left\| u _ {h} \right\| _ {0} ^ {2} + \left| u _ {h} \right| _ {1} ^ {2}, \tag {2.13}
$$

其中

$$
\left| u _ {h} \right| _ {1} ^ {2} = \sum_ {i = 1} ^ {N} h \left(\frac {u _ {i} - u _ {i - 1}}{h}\right) ^ {2}. \tag {2.14}
$$

若不特别说明, 我们用 $\|\cdot\|$ 表示 (2.11)-(2.13) 中任一种范数.

很明显, 要想差分解 $u_{h}$ 按范数 $\| \cdot \|$ 收敛到 $u$ , 差分算子 $L_{h}$ 必须在某种意义下逼近微分算子 $L$ , 这导致下列定义.

定义2.1 设 $\mathcal{M}$ 是某一充分光滑的函数类, $R_{h}(u)$ 是由截断误差(2.8)定义的网格函数. 若对任何 $u\in \mathcal{M}$ , 恒有

$$
\lim  _ {h \rightarrow 0} \| R _ {h} (u) \| = 0, \tag {2.15}
$$

则说差分算子 $L_{h}$ 逼近微分算子 $L$ ，而称(2.15）为相容条件

由（2.5）便知，差分算子（2.6）逼近微分算子(2.1)，且逼近的阶是 $\| R_h(u)\| _C = O\left(h^2\right),\| R_h(u)\| _0 = O\left(h^2\right),\| R_h(u)\| _1 = O(h).$

定义2.2 称差分解 $u_{h}$ 收敛到边值问题的解 $u$ ，如果当 $h$ 充分小时，(2.9)(2.10)的

解 $u_{h}$ 存在，且按某一范数 $\| \cdot \|$ 有

$$
\lim  _ {h \rightarrow 0} \| u _ {h} - u \| = 0, \tag {2.16}
$$

这里把 $u$ 看成 $I_h$ 上的网格函数.

将(2.4)写成

$$
L _ {h} u \left(x _ {i}\right) = f _ {i} + R _ {i} (u),
$$

以此与 (2.9) 相减, 得

$$
L _ {h} \left(u \left(x _ {i}\right) - u _ {i}\right) = R _ {i} (u).
$$

引进误差

$$
e _ {i} = u \left(x _ {i}\right) - u _ {i},
$$

则误差函数 $e_h(x_i) = e_i$ 满足下列差分方程：

$$
\left\{ \begin{array}{l} L _ {h} e _ {i} = R _ {i} (u), \\ e _ {0} = e _ {N} = 0, \end{array} \quad i = 1, 2, \dots , N - 1. \right. \tag {2.17}
$$

于是收敛性及收敛速度的估计问题，就归结到通过右端 $R_{i}(u)$ （截断误差）估计误差函数 $e_h$ 的问题.这和差分方程的稳定性有关.

定义2.3 称差分方程 $L_{h}v_{i} = f_{i}(i = 1,2,\dots ,N - 1),v_{0} = v_{N} = 0$ 关于右端稳定，如果存在与网格 $I_{h}$ 及右端 $f_{h}\left(f_{h}\left(x_{i}\right) = f_{i}\right)$ 无关的正常数 $M$ 和 $h_0$ ，使

$$
\left\| v _ {h} \right\| \leqslant M \left\| f _ {h} \right\| _ {R}, \quad 0 <   h <   h _ {0}, \tag {2.18}
$$

其中 $\| f_h\| _R$ 是右端 $f_{h}$ 的某一范数，它可以和 $\parallel \cdot \parallel$ 相同，也可以不同， $v_{h}(x_{i}) = v_{i},$ （20 $i = 1,2,\dots ,N - 1.$

不等式 (2.18) 表明, 解 $v_{h}$ 连续依赖右端 $f_{h}$ , 即右端变化小时解的变化也小. 实际上, 设 $u_{h}^{(1)}, u_{h}^{(2)}$ 是差分方程 (2.9) (2.10) 对应右端 $f_{h}^{(1)}, f_{h}^{(2)}$ 的解, 则 $v_{h} = u_{h}^{(1)} - u_{h}^{(2)}$ 满足 $L_{h} v_{i} = f_{i}^{(1)} - f_{i}^{(2)}, v_{0} = v_{N} = 0$ . 由 (2.18),

$$
\left\| u _ {h} ^ {(1)} - u _ {h} ^ {(2)} \right\| \leqslant M \left\| f _ {h} ^ {(1)} - f _ {h} ^ {(2)} \right\| _ {R}.
$$

若 $\lim_{h\to 0}\left\| f_h^{(1)} - f_h^{(2)}\right\| _R = 0$ ，则 $\lim_{h\to 0}\left\| u_h^{(1)} - u_h^{(2)}\right\| _R = 0.$

由（2.18）推出，与（2.9）（2.10）相应的齐方程 $(f_{i} = 0,\alpha = \beta = 0)$ 只能有平凡解 $u_{i}\equiv 0(i = 1,2,\dots ,N - 1),$ ，从而非齐方程对任何边值及右端有唯一解

将不等式 (2.18) 用到误差方程 (2.17), 则

$$
\left\| e _ {h} \right\| \leqslant M \left\| R _ {h} (u) \right\| _ {R}. \tag {2.19}
$$

若解 $u$ 充分光滑, $L_{h}$ 关于范数 $\| \cdot \| _R$ 满足相容条件, 则当 $h\to 0$ 时 $\| e_h\| \to 0$ ，从而差分

解收敛到边值问题的解，且有和截断误差相同的收敛阶

定理2.1 若边值问题的解 $u$ 充分光滑，差分方程按 $\| \cdot \| _R$ 满足相容条件，且关于右端稳定，则差分解 $u_{h}$ 按 $\| \cdot \|$ 收敛到边值问题的解，且有和 $\| R_h(u)\| _R$ 相同的收敛阶.

这样，为了建立差分解的收敛性，就需要检验相容条件并建立差分方程的稳定性。检验相容条件并不困难，例如由（2.9）定义的差分算子，我们曾用Taylor展式证明它关于 $\| \cdot \| _0$ 及 $\| \cdot \| _C$ 都满足相容条件，并且估计了截断误差的阶。我们的主要问题是去建立差分方程的稳定性，即建立形如（2.18）的估计式，称之为关于差分方程解的先验估计。

稳定性概念在理论研究和实际应用中都有重要意义. 实际上, 由于有实测误差和舍入误差, 右端数据不可能准确给出. 如果小的右端误差会引起解的很大偏离, 即差分方程不稳定, 便失去实际意义.

# 2.2 一维差分格式

考虑两点边值问题：

$$
L u = - \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) + r \frac {\mathrm {d} u}{\mathrm {d} x} + q u = f, \quad a <   x <   b, \tag {2.20}
$$

$$
u (a) = \alpha , \quad u (b) = \beta . \tag {2.21}
$$

假定 $p \in C^{1}[a, b], p(x) \geqslant p_{\min} > 0, r, q, f \in C[a, b], \alpha, \beta$ 是给定的常数.

本节我们将介绍构造差分格式的两种方法：直接差分化法、有限体积法. 还将讨论边值条件的逼近方法.

显然，我们可以造出许多逼近（2.20）（2.21）的差分格式，但并非任何格式都是可取的。一个好的差分格式，应该是以尽可能小的工作量（包括程序的准备和计算机的运算）获得所需精度的结果。因此，一方面，差分格式应该结构简单、便于求解；另一方面，应具有尽可能高的逼近阶。此外，还要根据问题的特点，对差分格式提出其他要求。

# 2.2.1 直接差分化

首先取 $N + 1$ 个节点：

$$
a = x _ {0} <   x _ {1} <   \dots <   x _ {i} <   \dots <   x _ {N} = b,
$$

将区间 $I = [a,b]$ 分成 $N$ 个小区间（见图2.1）：

$$
I _ {i}: x _ {i - 1} \leqslant x \leqslant x _ {i}, \quad i = 1, 2, \dots , N.
$$

于是得到 $I$ 的一个网格剖分. 记 $h_i = x_i - x_{i-1}$ , $h = \max_i h_i$ 为最大网格步长. 用 $I_h$ 表示网格内点 $x_1, x_2, \dots, x_{N-1}$ 的集合, $\bar{I}_h$ 表示内点和界点 $x_0 = a$ , $x_N = b$ 的集合.

![](images/1732bc0090debcfff62f18a17457eef8a6d32447244299e96343161cee87066f.jpg)  
图2.1 网格剖分

取相邻节点 $x_{i - 1},x_i$ 的中点 $x_{i - \frac{1}{2}} = \frac{1}{2} (x_{i - 1} + x_i)(i = 1,2,\dots ,N)$ ，称其为半整数点.则由节点

$$
a = x _ {0} <   x _ {\frac {1}{2}} <   x _ {\frac {3}{2}} <   \dots <   x _ {i - \frac {1}{2}} <   \dots <   x _ {N - \frac {1}{2}} <   x _ {N} = b
$$

又构成 $[a, b]$ 的一个剖分，称为对偶剖分。图2.1中打“·”号的是原剖分节点，打“×”号的是对偶剖分节点。

其次用差商代替微商, 将方程 (2.1) 在内点 $x_{i}$ 离散化. 注意对充分光滑的 $u$ , 由 Taylor 展开式有

$$
\frac {u \left(x _ {i + 1}\right) - u \left(x _ {i - 1}\right)}{h _ {i} + h _ {i + 1}} = \left[ \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i} + \frac {h _ {i + 1} - h _ {i}}{2} \left[ \frac {\mathrm {d} ^ {2} u}{\mathrm {d} x ^ {2}} \right] _ {i} + O \left(h ^ {2}\right), \tag {2.22}
$$

$$
\begin{array}{l} p \left(x _ {i - \frac {1}{2}}\right) \frac {u (x _ {i}) - u (x _ {i - 1})}{h _ {i}} = \left[ p \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i - \frac {1}{2}} + \frac {h _ {i} ^ {2}}{2 4} \left[ p \frac {\mathrm {d} ^ {3} u}{\mathrm {d} x ^ {3}} \right] _ {i - \frac {1}{2}} + O (h ^ {3}) \\ = \left[ p \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i - \frac {1}{2}} + \frac {h _ {i} ^ {2}}{2 4} \left[ p \frac {\mathrm {d} ^ {3} u}{\mathrm {d} x ^ {3}} \right] _ {i} + O (h ^ {3}), \tag {2.23} \\ \end{array}
$$

$$
p \left(x _ {i + \frac {1}{2}}\right) \frac {u \left(x _ {i + 1}\right) - u \left(x _ {i}\right)}{h _ {i + 1}} = \left[ p \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i + \frac {1}{2}} + \frac {h _ {i + 1} ^ {2}}{2 4} \left[ p \frac {\mathrm {d} ^ {3} u}{\mathrm {d} x ^ {3}} \right] _ {i} + O \left(h ^ {3}\right). \tag {2.24}
$$

由(2.24)减(2.23)，并除以 $\frac{(h_i + h_{i + 1})}{2}$ 得

$$
\begin{array}{l} \frac {2}{h _ {i} + h _ {i + 1}} \left[ p \left(x _ {i + \frac {1}{2}}\right) \frac {u \left(x _ {i + 1}\right) - u \left(x _ {i}\right)}{h _ {i + 1}} - p \left(x _ {i - \frac {1}{2}}\right) \frac {u \left(x _ {i}\right) - u \left(x _ {i - 1}\right)}{h _ {i}} \right] \\ = \frac {2}{h _ {i} + h _ {i + 1}} \left(\left[ p \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i + \frac {1}{2}} - \left[ p \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i - \frac {1}{2}}\right) + \frac {h _ {i + 1} - h _ {i}}{1 2} \left[ p \frac {\mathrm {d} ^ {3} u}{\mathrm {d} x ^ {3}} \right] _ {i} + O (h ^ {2}) \\ = \left[ \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) \right] _ {i} + \frac {h _ {i + 1} - h _ {i}}{4} \left[ \frac {\mathrm {d} ^ {2}}{\mathrm {d} x ^ {2}} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) \right] _ {i} + \frac {h _ {i + 1} - h _ {i}}{1 2} \left[ p \frac {\mathrm {d} ^ {3} u}{\mathrm {d} x ^ {3}} \right] _ {i} + O (h ^ {2}). \tag {2.25} \\ \end{array}
$$

令 $p_{i - \frac{1}{2}} = p\left(x_{i - \frac{1}{2}}\right),r_i = r\left(x_i\right),q_i = q\left(x_i\right),f_i = f\left(x_i\right)$ ，则由（2.22）（2.25）知，边值问题的解 $u(x)$ 满足方程：

$$
\begin{array}{l} L _ {h} u \left(x _ {i}\right) \equiv - \frac {2}{h _ {i} + h _ {i + 1}} \left[ p _ {i + \frac {1}{2}} \frac {u \left(x _ {i + 1}\right) - u \left(x _ {i}\right)}{h _ {i + 1}} - p _ {i - \frac {1}{2}} \frac {u \left(x _ {i}\right) - u \left(x _ {i - 1}\right)}{h _ {i}} \right] + \\ \frac {r _ {i}}{h _ {i} + h _ {i + 1}} \left[ u \left(x _ {i + 1}\right) - u \left(x _ {i - 1}\right) \right] + q _ {i} u \left(x _ {i}\right) \\ \end{array}
$$

$$
= f _ {i} + R _ {i} (u), \tag {2.26}
$$

其中

$$
R _ {i} (u) = - \left(h _ {i + 1} - h _ {i}\right) \left(\frac {1}{4} \left[ \frac {\mathrm {d} ^ {2}}{\mathrm {d} x ^ {2}} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) \right] _ {i} + \frac {1}{1 2} \left[ p \frac {\mathrm {d} ^ {3} u}{\mathrm {d} x ^ {3}} \right] _ {i} - \frac {1}{2} \left[ r \frac {\mathrm {d} ^ {2} u}{\mathrm {d} x ^ {2}} \right] _ {i}\right) + O \left(h ^ {2}\right)
$$

为差分算子 $L_{h}$ 的截断误差. 舍去 $R_{i}(u)$ , 便得逼近边值问题 (2.20) (2.21) 的差分方程:

$$
L _ {h} u _ {i} \equiv - \frac {2}{h _ {i} + h _ {i + 1}} \left[ p _ {i + \frac {1}{2}} \frac {u _ {i + 1} - u _ {i}}{h _ {i + 1}} - p _ {i - \frac {1}{2}} \frac {u _ {i} - u _ {i - 1}}{h _ {i}} \right] +
$$

$$
\frac {r _ {i}}{h _ {i} + h _ {i + 1}} \left(u _ {i + 1} - u _ {i - 1}\right) + q _ {i} u _ {i} = f _ {i}, \quad i = 1, 2, \dots , N - 1, \tag {2.28}
$$

$$
u _ {0} = \alpha , \quad u _ {N} = \beta . \tag {2.29}
$$

差分方程 (2.28) 也可用数值微商公式

$$
\begin{array}{l} \left[ \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i} \approx \frac {u _ {i + 1} - u _ {i - 1}}{h _ {i} + h _ {i + 1}}, \\ \left[ \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) \right] _ {i} \approx \left(p _ {i + \frac {1}{2}} \left[ \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i + \frac {1}{2}} - p _ {i - \frac {1}{2}} \left[ \frac {\mathrm {d} u}{\mathrm {d} x} \right] _ {i - \frac {1}{2}}\right) / \frac {h _ {i} + h _ {i + 1}}{2} \\ \approx \frac {2}{h _ {i} + h _ {i + 1}} \left(p _ {i + \frac {1}{2}} \frac {u _ {i + 1} - u _ {i}}{h _ {i + 1}} - p _ {i - \frac {1}{2}} \frac {u _ {i} - u _ {i - 1}}{h _ {i}}\right) \\ \end{array}
$$

代入方程 (2.20) 得到. 我们采取前述推导, 是为了导出截断误差公式 (2.27).

方程 (2.28) (2.29) 是 $N - 1$ 阶的线性代数方程组。若节点次序由左到右排列，则系数矩阵 $\mathbf{A}$ 是三对角矩阵。由于 $r$ 不恒等于零，矩阵 $\mathbf{A}$ 不对称。当 $r \equiv 0$ 即 (2.20) 对称时，若网格不均匀，则矩阵 $\mathbf{A}$ 也可能不对称，但可以对称化，这只要在 (2.28) 两端乘 $(h_i + h_{i+1})$ 即可。求解 (2.28) (2.29) 就得出解 $u(x)$ 在 $x_i$ 的近似值 $u_i$ 。

由方程 (2.26) (2.28), 截断误差 $R_{i}(u)$ 可表为

$$
R _ {i} (u) = L _ {h} u \left(x _ {i}\right) - L _ {h} u _ {i} = L _ {h} \left(u \left(x _ {i}\right) - u _ {i}\right). \tag {2.30}
$$

以 $R_{h}(u)$ 表示由 $R_{i}(u)$ 定义的网格函数, 则由 (2.27) 可知截断误差按 $\| \cdot \|_{C}$ 或 $\| \cdot \|_{0}$ 的阶都是 $O(h)$ . 当网格均匀, 即 $h_i = h(i = 1,2,\dots ,N)$ 时, $\| R_h(u)\| _c$ 或 $\| R_h(u)\| _0$ 的阶提高为 $O\left(h^{2}\right)$ . 此时差分方程 (2.28) 简化为

$$
\begin{array}{l} L _ {h} u _ {i} = - \frac {1}{h ^ {2}} \left[ p _ {i + \frac {1}{2}} u _ {i + 1} - \left(p _ {i + \frac {1}{2}} + p _ {i - \frac {1}{2}}\right) u _ {i} + p _ {i - \frac {1}{2}} u _ {i - 1} \right] + \\ r _ {i} \frac {u _ {i + 1} - u _ {i - 1}}{2 h} + q _ {i} u _ {i} = f _ {i}. \tag {2.31} \\ \end{array}
$$

这相当于用一阶中心差商、二阶中心差商依次代替方程(2.20)的一阶微商和二阶微商

# 2.2.2 有限体积法（积分插值法）

考虑守恒型微分方程：

$$
L u = - \frac {\mathrm {d}}{\mathrm {d} x} \left(p (x) \frac {\mathrm {d} u}{\mathrm {d} x}\right) + q (x) u = f (x). \tag {2.32}
$$

如果把它看作是分布在一根杆上的稳定温度场方程，则在 $[a,b]$ 上任一小区间 $\left[x^{(1)},x^{(2)}\right]$ 上的热量守恒律具有形式

$$
- \int_ {x ^ {(1)}} ^ {x ^ {(2)}} \frac {\mathrm {d}}{\mathrm {d} x} \left(p (x) \frac {\mathrm {d} u}{\mathrm {d} x}\right) \mathrm {d} x + \int_ {x ^ {(1)}} ^ {x ^ {(2)}} q u \mathrm {d} x = \int_ {x ^ {(1)}} ^ {x ^ {(2)}} f \mathrm {d} x,
$$

或

$$
W \left(x ^ {(1)}\right) - W \left(x ^ {(2)}\right) + \int_ {x ^ {(1)}} ^ {x ^ {(2)}} q (x) u \mathrm {d} x = \int_ {x ^ {(1)}} ^ {x ^ {(2)}} f \mathrm {d} x, \tag {2.33}
$$

其中

$$
W (x) = p (x) \frac {\mathrm {d} u}{\mathrm {d} x}. \tag {2.34}
$$

把微分方程 (2.32) 写成积分守恒型 (2.33) 后, 最高阶微商由二阶降到一阶, 从而可减弱对 $p$ 、 $u$ 光滑性的要求. 以后会看到, 从积分守恒型方程出发构造差分格式, 便于推广到任意网格和处理第二边值条件.

既然具守恒形式的微分方程反映了物理、力学某些守恒定律，那么我们构造的差分格式也应反映这一基本性质.现在来构造这种差分格式

特别地, 于 (2.33) 中取 $\left[x^{(1)}, x^{(2)}\right]$ 为对偶单元 $\left[x_{i - \frac{1}{2}}, x_{i + \frac{1}{2}}\right]$ , 则

$$
W \left(x _ {i - \frac {1}{2}}\right) - W \left(x _ {i + \frac {1}{2}}\right) + \int_ {x _ {i - \frac {1}{2}}} ^ {x _ {i + \frac {1}{2}}} q u \mathrm {d} x = \int_ {x _ {i - \frac {1}{2}}} ^ {x _ {i + \frac {1}{2}}} f \mathrm {d} x. \tag {2.35}
$$

考虑到 $p(x)$ 可能有间断点，此时由(2.34)进一步差分化是不合适的.但“热流量” $W(x)$ 恒连续，故将(2.34)改写成

$$
\frac {\mathrm {d} u}{\mathrm {d} x} = \frac {W (x)}{p (x)},
$$

再沿 $[x_{i - 1},x_i]$ 积分，得

$$
u _ {i} - u _ {i - 1} = \int_ {x _ {i - 1}} ^ {x _ {i}} \frac {W (x)}{p (x)} \mathrm {d} x,
$$

利用中矩形公式，有

$$
W _ {i - \frac {1}{2}} \approx a _ {i} \frac {u _ {i} - u _ {i - 1}}{h _ {i}}, \tag {2.36}
$$

$$
a _ {i} = \left(\frac {1}{h _ {i}} \int_ {x _ {i - 1}} ^ {x _ {i}} \frac {\mathrm {d} x}{p (x)}\right) ^ {- 1}. \tag {2.37}
$$

又

$$
\int_ {x _ {i - \frac {1}{2}}} ^ {x _ {i + \frac {1}{2}}} q u \mathrm {d} x \approx \frac {h _ {i} + h _ {i + 1}}{2} d _ {i} u _ {i}, \tag {2.38}
$$

$$
d _ {i} = \frac {2}{h _ {i} + h _ {i + 1}} \int_ {x _ {i - \frac {1}{2}}} ^ {x _ {i + \frac {1}{2}}} q (x) \mathrm {d} x. \tag {2.39}
$$

将(2.36)(2.38)代到(2.35)，即得守恒型差分方程：

$$
- \left(a _ {i + 1} \frac {u _ {i + 1} - u _ {i}}{h _ {i + 1}} - a _ {i} \frac {u _ {i} - u _ {i - 1}}{h _ {i}}\right) + \frac {1}{2} \left(h _ {i} + h _ {i + 1}\right) d _ {i} u _ {i} = \frac {1}{2} \left(h _ {i} + h _ {i + 1}\right) \phi_ {i}, \tag {2.40}
$$

$$
\phi_ {i} = \frac {2}{h _ {i} + h _ {i + 1}} \int_ {x _ {i - \frac {1}{2}}} ^ {x _ {i + \frac {1}{2}}} f (x) d x. \tag {2.41}
$$

如果系数 $p, q$ 及右端 $f$ 光滑，则可用中矩形公式计算 (2.37) (2.39) 和 (2.41)，从而

$$
\left\{ \begin{array}{l} a _ {i} = p _ {i - \frac {1}{2}} = p \left(x _ {i - \frac {1}{2}}\right), \\ d _ {i} = q _ {i} = q \left(x _ {i}\right), \\ \phi_ {i} = f _ {i} = f \left(x _ {i}\right), \end{array} \right. \tag {2.42}
$$

也可用梯形公式，此时

$$
\left\{ \begin{array}{l} a _ {i} = \frac {2 p _ {i - 1} p _ {i}}{p _ {i - 1} + p _ {i}}, \\ d _ {i} = \frac {1}{2} \left(q _ {i - \frac {1}{2}} + q _ {i + \frac {1}{2}}\right), \\ f _ {i} = \frac {1}{2} \left(f _ {i - \frac {1}{2}} + f _ {i + \frac {1}{2}}\right). \end{array} \right. \tag {2.43}
$$

注2.1差分方程（2.40）也适用于具第一类间断系数的微分方程.此时若系数计算公式(2.42)或(2.43)右端的函数在间断点取值，则应取左右极限的算术平均值.对于具间断系数的微分方程，保持守恒形式尤为重要.例如微分方程

$$
\frac {\mathrm {d}}{\mathrm {d} x} \left(p (x) \frac {\mathrm {d} u}{\mathrm {d} x}\right) = 0.
$$

若把它写成非守恒形式

$$
p (x) \frac {\mathrm {d} ^ {2} u}{\mathrm {d} x ^ {2}} + p ^ {\prime} (x) \frac {\mathrm {d} u}{\mathrm {d} x} = 0,
$$

再用中心差分格式

$$
p _ {i} \frac {u _ {i + 1} - 2 u _ {i} + u _ {i - 1}}{h ^ {2}} + \frac {p _ {i + 1} - p _ {i - 1}}{2 h} \frac {u _ {i + 1} - u _ {i - 1}}{2 h} = 0,
$$

则差分解可能不收敛 (参看 [16] 的第 2 章 §2).

# 2.2.3 边值条件的处理

最简单的第一边值条件已处理过了，现在处理第二、第三边值条件：

$$
u ^ {\prime} (a) = \alpha_ {0} u (a) + \alpha_ {1}, \tag {2.44}
$$

$$
u ^ {\prime} (b) = \beta_ {0} u (b) + \beta_ {1}. \tag {2.45}
$$

最容易想到的是用数值微商公式

$$
u ^ {\prime} (a) \approx \frac {u _ {1} - u _ {0}}{h _ {1}}, \quad u ^ {\prime} (b) \approx \frac {u _ {N} - u _ {N - 1}}{h _ {N}}
$$

代替 (2.44) 和 (2.45) 中的微商. 但这样处理有两个缺点: 一是截断误差的阶比内点低, 例如对均匀网格, 内点的截断误差为 $O(h^{2})$ , 界点的截断误差的阶为 $O(h)$ ; 二是可能会破坏差分方程 (2.40) 的对称性. 为此我们用有限体积法, 像推导内点差分方程那样导出近似边值条件.

因为 $p(x) > 0$ ，不失一般性可将边值条件(2.44）和（2.45）写成形式

$$
- p (a) u ^ {\prime} (a) = \alpha_ {0} u (a) + \alpha_ {1}, \tag {2.46}
$$

$$
- p (b) u ^ {\prime} (b) = \beta_ {0} u (b) + \beta_ {1}. \tag {2.47}
$$

于积分守恒形式 (2.33) 中取 $x^{(1)} = x_0 = a, x^{(2)} = x_{\frac{1}{2}}$ 得

$$
W (a) - W \left(x _ {\frac {1}{2}}\right) + \int_ {x _ {0}} ^ {x _ {\frac {1}{2}}} q u d x = \int_ {x _ {0}} ^ {x _ {\frac {1}{2}}} f d x.
$$

而

$$
W (a) = \left. p (x) \frac {\mathrm {d} u}{\mathrm {d} x} \right| _ {x = a} = - \left(\alpha_ {0} u _ {0} + \alpha_ {1}\right),
$$

故

$$
- W \left(x _ {\frac {1}{2}}\right) + \int_ {x _ {0}} ^ {x _ {\frac {1}{2}}} q u d x = \left(\alpha_ {0} u _ {0} + \alpha_ {1}\right) + \int_ {x _ {0}} ^ {x _ {\frac {1}{2}}} f d x. \tag {2.48}
$$

由(2.36)得

$$
W \left(x _ {\frac {1}{2}}\right) \approx a _ {1} \frac {u _ {1} - u _ {0}}{h _ {1}}, a _ {1} = \left(\frac {1}{h _ {1}} \int_ {x _ {0}} ^ {x _ {1}} \frac {\mathrm {d} x}{p (x)}\right) ^ {- 1}. \tag {2.49}
$$

又

$$
\int_ {x _ {0}} ^ {x _ {\frac {1}{2}}} q u \mathrm {d} x \approx \frac {h _ {1}}{2} d _ {0} u _ {0}, d _ {0} = \frac {2}{h _ {1}} \int_ {x _ {0}} ^ {x _ {\frac {1}{2}}} q \mathrm {d} x, \tag {2.50}
$$

$$
\int_ {x _ {0}} ^ {x \frac {1}{2}} f \mathrm {d} x = \frac {h _ {1}}{2} \phi_ {0}, \phi_ {0} = \frac {2}{h _ {1}} \int_ {x _ {0}} ^ {x \frac {1}{2}} f \mathrm {d} x. \tag {2.51}
$$

以(2.49)一(2.51)代到(2.48)，得

$$
- a _ {1} \frac {u _ {1} - u _ {0}}{h _ {1}} + \left(- \alpha_ {0} + \frac {h _ {1}}{2} d _ {0}\right) u _ {0} - \left(\alpha_ {1} + \frac {h _ {1}}{2} \phi_ {0}\right) = 0. \tag {2.52}
$$

类似地可导出逼近（2.47）的差分方程

可以证明, 当网格均匀且系数光滑时, 差分方程 (2.40) 逼近 (2.32) 的阶为 $O(h^{2})$ , 边界差分方程 (2.52) 逼近 (2.44) 的阶亦为 $O(h^{2})$

# 2.2.4 习题

1. 用有限体积法导出逼近微分方程 (2.20) 的差分方程  
2.构造逼近

$$
(p u ^ {\prime \prime}) ^ {\prime \prime} + (q u ^ {\prime}) ^ {\prime} + r u = f,   \text {于} (a, b),
$$

$$
u (a) = u ^ {\prime} (a) = 0, \quad u (b) = u ^ {\prime} (b) = 0
$$

的中心差分格式

# 2.3 矩形网的差分格式

考虑 Poisson 方程:

$$
- \Delta u = f (x, y), \quad (x, y) \in G, \tag {2.53}
$$

其中 $G$ 是 $xy$ 平面上一有界区域, 其边界 $\Gamma$ 为分段光滑曲线. 在 $\Gamma$ 上 $u$ 满足下列边值条件之一:

$$
u | _ {\Gamma} = \alpha (x, y) \quad (\text {第 一 边 值 条 件}), \tag {2.54}
$$

$$
\left. \frac {\partial u}{\partial n} \right| _ {\Gamma} = \beta (x, y) \quad (\text {第 二 边 值 条 件}), \tag {2.55}
$$

$$
\frac {\partial u}{\partial \boldsymbol {n}} + k u | _ {\Gamma} = \gamma (\boldsymbol {x}, y) \quad (\text {第 三 边 值 条 件}), \tag {2.56}
$$

其中 $f(x,y),\alpha (x,y),\beta (x,y),\gamma (x,y)$ 及 $k(x,y)$ 都是连续函数， $k\geqslant 0$ .本节讨论逼近方程(2.53）及相应边值条件的差分格式.原则上，前节的方法都可推广到二维边值问题，但会遇到新的困难.例如随着维数的增加，求解域的几何形状会更复杂，如何作网格剖分及处理好边值条件，就是一个重要而困难的问题

# 2.3.1 五点差分格式

取定沿 $x$ 轴和 $y$ 轴的步长 $h_1$ 和 $h_2$ , $h = \left(h_1^2 + h_2^2\right)^{\frac{1}{2}}$ . 作两族与坐标轴平行的直线:

$$
x = i h _ {1}, \quad i = 0, \pm 1, \dots ,
$$

$$
y = j h _ {2}, \quad j = 0, \pm 1, \dots .
$$

两族直线的交点 $(ih_{1},jh_{2})$ 称为网点或节点，记为 $(x_{i},y_{j})$ 或 $(i,j)$ 两个节点 $(x_{i},y_{j})$ 和 $(x_{i^{\prime}},y_{j^{\prime}})$ 说是相邻，如果

$$
\left| \frac {x _ {i} - x _ {i ^ {\prime}}}{h _ {1}} \right| + \left| \frac {y _ {j} - y _ {j ^ {\prime}}}{h _ {2}} \right| = 1 \quad \text {或} \quad | i - i ^ {\prime} | + | j - j ^ {\prime} | = 1.
$$

以 $G_{h} = \{(x_{i},y_{i})\in G\}$ 表示所有属于 $G$ 内部的节点集合，并称如此的节点为内点. 以 $\Gamma_{h}$ 表示网线 $x = x_{i}$ 或 $y = y_{j}$ 与 $\Gamma$ 的交点集合，并称如此的点为界点. 令 $\bar{G}_h = G_h\cup \Gamma_h$ 则 $\bar{G}_h$ 就是代替域 $\bar{G} = G\cup \Gamma$ 的网点集合. 若内点 $(x_{i},y_{j})$ 的四个相邻点都属于 $G_{h}$ ，就称为正则内点；否则称为非正则内点. 图2.2中打“○”号的点为正则内点，打“×”号的点为非正则内点，打“·”号的点为界点.

![](images/b51c81c4f96f97dea131f1261d04a87a960b74fdc781721b2ecca7e30f553b6b.jpg)  
图2.2 不同节点类型

现在假定 $(x_{i},y_{j})$ 为正则内点.沿 $x,y$ 方向分别用二阶中心差商代替 $u_{xx},u_{yy}$ ，则得

$$
- \Delta_ {h} u _ {i j} = - \left(\frac {u _ {i + 1 , j} - 2 u _ {i j} + u _ {i - 1 , j}}{h _ {1} ^ {2}} + \frac {u _ {i , j + 1} - 2 u _ {i j} + u _ {i , j - 1}}{h _ {2} ^ {2}}\right) = f _ {i j}, \tag {2.57}
$$

式中 $u_{ij}$ 表示节点 $(i,j)$ 上的网函数值. 若以 $u_h, f_h$ 表示网格函数， $u_h(x_i, y_j) = u_{ij}, f_h(x_i, y_j) = f_{ij} = f(x_i, y_j)$ ，则差分方程 (2.57) 可简写成

$$
- \Delta_ {h} u _ {h} = f _ {h}. \tag {2.58}
$$

利用Taylor展开式，

$$
\begin{array}{l} \frac {u \left(x _ {i + 1} , y _ {j}\right) - 2 u \left(x _ {i} , y _ {j}\right) + u \left(x _ {i - 1} , y _ {j}\right)}{h _ {1} ^ {2}} \\ = \frac {\partial^ {2} u \left(x _ {i} , y _ {j}\right)}{\partial x ^ {2}} + \frac {h _ {1} ^ {2}}{1 2} \frac {\partial^ {4} u \left(x _ {i} , y _ {j}\right)}{\partial x ^ {4}} + \frac {h _ {1} ^ {4}}{3 6 0} \frac {\partial^ {6} u \left(x _ {i} , y _ {j}\right)}{\partial x ^ {6}} + O \left(h _ {1} ^ {6}\right), (2.59) \\ \frac {u \left(x _ {i} , y _ {j + 1}\right) - 2 u \left(x _ {i} , y _ {j}\right) + u \left(x _ {i} , y _ {j - 1}\right)}{h _ {2} ^ {2}} \\ = \frac {\partial^ {2} u (x _ {i} , y _ {j})}{\partial y ^ {2}} + \frac {h _ {2} ^ {2}}{1 2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial y ^ {4}} + \frac {h _ {2} ^ {4}}{3 6 0} \frac {\partial^ {6} u (x _ {i} , y _ {j})}{\partial y ^ {6}} + O \left(h _ {2} ^ {6}\right), (2.60) \\ \end{array}
$$

可得差分算子 $-\Delta_{h}$ 的截断误差

$$
\begin{array}{l} R _ {i j} (u) = \Delta u \left(x _ {i}, y _ {j}\right) - \Delta_ {h} u \left(x _ {i}, y _ {j}\right) \\ = - \frac {1}{1 2} \left[ h _ {1} ^ {2} \frac {\partial^ {4} u}{\partial x ^ {4}} + h _ {2} ^ {2} \frac {\partial^ {4} u}{\partial y ^ {4}} \right] _ {i j} + O (h ^ {4}) \\ = O \left(h ^ {2}\right), \tag {2.61} \\ \end{array}
$$

其中 $u$ 是方程(2.53)的光滑解

由于差分方程 (2.57) 中只出现 $u$ 在 $(i,j)$ 及其四个邻点上的值, 故称为五点差分格式, 其图式如图 2.3. 特别地, 取正方形网格: $h_1 = h_2 = h$ , 则差分方程 (2.57) 简化为

$$
u _ {i j} - \frac {1}{4} \left(u _ {i - 1, j} + u _ {i, j - 1} + u _ {i + 1, j} + u _ {i, j + 1}\right) = \frac {h ^ {2}}{4} f _ {i j}. \tag {2.62}
$$

若 $f\equiv 0$ (Laplace方程),则

$$
u _ {i j} = \frac {1}{4} \left(u _ {i - 1, j} + u _ {i, j - 1} + u _ {i + 1, j} + u _ {i, j + 1}\right). \tag {2.63}
$$

![](images/389f9a55308d85632b75432386a2c08640211b6fb6831e8a5e529ea362b951d7.jpg)  
图2.3 五点差分格式

注2.2 若将(2.59)(2.60)两式相加，则得

$$
\begin{array}{l} \Delta_ {h} u \left(x _ {i}, y _ {j}\right) \\ = \Delta u (x _ {i}, y _ {j}) + \frac {1}{1 2} \left(h _ {1} ^ {2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {4}} + h _ {2} ^ {2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial y ^ {4}}\right) + O (h ^ {4}) \\ \end{array}
$$

$$
\begin{array}{l} = \Delta u (x _ {i}, y _ {j}) + \frac {1}{1 2} \left(h _ {1} ^ {2} \frac {\partial^ {2}}{\partial x ^ {2}} + h _ {2} ^ {2} \frac {\partial^ {2}}{\partial y ^ {2}}\right) \left(\frac {\partial^ {2} u (x _ {i} , y _ {j})}{\partial x ^ {2}} + \frac {\partial^ {2} u (x _ {i} , y _ {j})}{\partial y ^ {2}}\right) - \\ \frac {h _ {1} ^ {2} + h _ {2} ^ {2}}{1 2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {2} \partial y ^ {2}} + O (h ^ {4}) \\ = - f \left(x _ {i}, y _ {j}\right) - \frac {1}{1 2} \left(h _ {1} ^ {2} \frac {\partial^ {2} f \left(x _ {i} , y _ {j}\right)}{\partial x ^ {2}} + h _ {2} ^ {2} \frac {\partial^ {2} f \left(x _ {i} , y _ {j}\right)}{\partial y ^ {2}}\right) - \\ \frac {h _ {1} ^ {2} + h _ {2} ^ {2}}{1 2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {2} \partial y ^ {2}} + O (h ^ {4}). \\ \end{array}
$$

又

$$
\begin{array}{l} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {2} \partial y ^ {2}} = \frac {u _ {x x} ^ {\prime \prime} (x _ {i} , y _ {j + 1}) - 2 u _ {x x} ^ {\prime \prime} (x _ {i} , y _ {j}) + u _ {x x} ^ {\prime \prime} (x _ {i} , y _ {j - 1})}{h _ {2} ^ {2}} + O \left(h _ {2} ^ {2}\right) \\ = \frac {1}{h _ {1} ^ {2} h _ {2} ^ {2}} [ u (x _ {i + 1}, y _ {j + 1}) - 2 u (x _ {i}, y _ {j + 1}) + \\ u \left(x _ {i - 1}, y _ {j + 1}\right) - 2 \left(u \left(x _ {i + 1}, y _ {j}\right) - 2 u \left(x _ {i}, y _ {j}\right) + \right. \\ u \left(x _ {i - 1}, y _ {j}\right)) + u \left(x _ {i + 1}, y _ {j - 1}\right) - \\ \left. 2 u \left(x _ {i}, y _ {j - 1}\right) + u \left(x _ {i - 1}, y _ {j - 1}\right) \right] + O \left(h ^ {2}\right). \\ \end{array}
$$

因此

$$
\begin{array}{l} \Delta_ {h} u \left(x _ {i}, y _ {j}\right) + \frac {1}{1 2} \left[ 4 u \left(x _ {i}, y _ {j}\right) - 2 \left(u \left(x _ {i - 1}, y _ {j}\right) + u \left(x _ {i}, y _ {j - 1}\right) + \right. \right. \\ u \left(x _ {i + 1}, y _ {j}\right) + u \left(x _ {i}, y _ {j + 1}\right)) + u \left(x _ {i - 1}, y _ {j - 1}\right) + u \left(x _ {i + 1}, y _ {j - 1}\right) + \\ \left. u \left(x _ {i + 1}, y _ {j + 1}\right) + u \left(x _ {i - 1}, y _ {j + 1}\right) \right] \frac {h _ {1} ^ {2} + h _ {2} ^ {2}}{h _ {1} ^ {2} h _ {2} ^ {2}} \\ = - f \left(x _ {i}, y _ {j}\right) - \frac {1}{1 2} \left(h _ {1} ^ {2} \frac {\partial^ {2} f \left(x _ {i} , y _ {j}\right)}{\partial x ^ {2}} + h _ {2} ^ {2} \frac {\partial^ {2} f \left(x _ {i} , y _ {j}\right)}{\partial y ^ {2}}\right) + O \left(h ^ {4}\right). \\ \end{array}
$$

舍去截断误差项，便得到逼近 Poisson 方程的九点差分格式：

$$
\begin{array}{l} - \Delta_ {h} u _ {i j} - \frac {1}{1 2} [ 4 u _ {i j} - 2 (u _ {i - 1, j} + u _ {i, j - 1} + u _ {i + 1, j} + u _ {i, j + 1}) + \\ \left. u _ {i - 1, j - 1} + u _ {i + 1, j - 1} + u _ {i + 1, j + 1} + u _ {i - 1, j + 1} \right] \frac {h _ {1} ^ {2} + h _ {2} ^ {2}}{h _ {1} ^ {2} h _ {2} ^ {2}} \\ = f _ {i j} + \frac {1}{1 2} \left[ h _ {1} ^ {2} f _ {x x} ^ {\prime \prime} (x _ {i}, y _ {j}) + h _ {2} ^ {2} f _ {y y} ^ {\prime \prime} (x _ {i}, y _ {j}) \right], \\ \end{array}
$$

其截断误差的阶为 $O\left(h^4\right)$

现在用有限体积法推导五点差分格式. 为此我们需要作对偶剖分. 记 $x_{i - \frac{1}{2}} = \left(i - \frac{1}{2}\right)h_1$ $y_{j - \frac{1}{2}} = \left(j - \frac{1}{2}\right)h_2.$ 作两族平行于坐标轴的直线 $x = x_{i - \frac{1}{2}}$ 和 $y = y_{j - \frac{1}{2}}, i, j = 0, \pm 1, \dots,$ 其交点属于 $G$ 内部者为对偶剖分的内点，直线与边界 $\Gamma$ 的交点为对偶剖分的界点. 对于任一正则内点 $(x_i, y_j)$ ，考虑对偶剖分的网点： $A\left(x_{i - \frac{1}{2}}, y_{j - \frac{1}{2}}\right), B\left(x_{i + \frac{1}{2}}, y_{j - \frac{1}{2}}\right)$ ，

$C\left(x_{i + \frac{1}{2}},y_{j + \frac{1}{2}}\right),D\left(x_{i - \frac{1}{2}},y_{j + \frac{1}{2}}\right)$ 用 $G_{ij}$ 表示以 $A,B,C,D$ 为顶点的矩形域，称为控制体积或对偶单元， $\partial G_{ij}$ 为其边界(参看图2.4).于 $G_{ij}$ 积分Poisson方程(2.53)，并利用Green公式，得到Poisson方程的积分守恒形式：

$$
- \int_ {\partial G _ {i j}} \frac {\partial u}{\partial \boldsymbol {n}} \mathrm {d} s = \iint_ {G _ {i j}} f \mathrm {d} x \mathrm {d} y. \tag {2.64}
$$

![](images/261fb0b6bdd44bcf5aa3cd71f099f7be5645fdde168af2f7ca870ce77df411be.jpg)  
图2.4 对偶剖分

式中 $\frac{\partial u}{\partial n}$ 表示 $u$ 沿矩形 $\partial G_{ij}$ 的外法向导数. 用中矩形公式代替沿四边的线积分, 再用中心差商代替外法向导数, 则

$$
\int_ {\partial G _ {i j}} \frac {\partial u}{\partial \boldsymbol {n}} d s \approx \frac {u _ {i , j - 1} - u _ {i j}}{h _ {2}} h _ {1} + \frac {u _ {i + 1 , j} - u _ {i j}}{h _ {1}} h _ {2} + \frac {u _ {i , j + 1} - u _ {i j}}{h _ {2}} h _ {1} + \frac {u _ {i - 1 , j} - u _ {i j}}{h _ {1}} h _ {2}.
$$

以之代到 (2.64), 并除以 $h_1 h_2$ , 即得五点差分格式:

$$
- \left(\frac {u _ {i + 1 , j} - 2 u _ {i j} + u _ {i - 1 , j}}{h _ {1} ^ {2}} + \frac {u _ {i , j + 1} - 2 u _ {i j} + u _ {i , j - 1}}{h _ {2} ^ {2}}\right) = \phi_ {i j}, \tag {2.65}
$$

$$
\phi_ {i j} = \frac {1}{h _ {1} h _ {2}} \iint_ {G _ {i j}} f d x d y \approx f _ {i j},
$$

它和 (2.57) 一致

# 2.3.2 边值条件的处理

先讨论第一边值条件

$$
u | _ {\Gamma} = \alpha (x, y). \tag {2.66}
$$

以 $G_h^*$ 表示非正则内点集合, $\Gamma_h$ 表示界点集合. 因 $\Gamma_h \subset \Gamma$ , 当 $(\bar{x}_i, \bar{y}_j) \in \Gamma_h$ 时, 便令

$$
u _ {i j} = \alpha (\bar {x} _ {i}, \bar {y} _ {j}). \tag {2.67}
$$

当 $(\bar{x}_i,\bar{y}_j)\in G_h^*$ 即它为非正则内点时，将它和正则内点一样，在 $(\bar{x}_i,\bar{y}_j)$ 建立逼近Poisson方程的不等距差分格式.例如在节点0(参看图2.5)有

$$
\left. - \left[ \frac {1}{\bar {h} _ {1}} \left(\frac {u _ {1} - u _ {0}}{h _ {1}} - \frac {u _ {0} - u _ {3}}{h _ {1} ^ {-}}\right) + \frac {1}{h _ {2} ^ {2}} \left(u _ {2} - 2 u _ {0} + u _ {4}\right) \right] = f _ {0}, \right. \tag {2.68}
$$

$$
\bar {h} _ {1} = \frac {1}{2} \left(h _ {1} + h _ {1} ^ {-}\right), h _ {1} ^ {-} = x _ {0} - x _ {3}.
$$

其截断误差的阶为 $O(h)$

![](images/40fe88d7a08a3c2dd3d5e918b607655099e86b4139491ae87e7f804fc19afbc6.jpg)  
图2.5 非正则内点

按(2.68)处理边值条件有一个缺点，即破坏了对称正定性，而这一性质是五点差分格式所固有的。为了保持对称正定性，可用

$$
- \left[ \frac {1}{h _ {1}} \left(\frac {u _ {1} - u _ {0}}{h _ {1}} - \frac {u _ {0} - u _ {3}}{h _ {1} ^ {-}}\right) + \frac {1}{h _ {2} ^ {2}} \left(u _ {2} - 2 u _ {0} + u _ {4}\right) \right] = f _ {0} \tag {2.69}
$$

代替 (2.68). 此时截断误差的阶降为 $O(1)$ . 尽管如此, 仍可证明, 差分解的收敛阶仍是 $O(h^{2})$ (按最大模) (见 [16] 的第 3 章).

今讨论第二、第三边值条件：

$$
\left. \frac {\partial u}{\partial n} + k u \right| _ {\Gamma} = \gamma . \tag {2.70}
$$

我们用有限体积法构造 (2.70) 的差分逼近. 我们只讨论一种特殊情况, 假定 $\Gamma_h$ 中的节点 (界点) 是两族网线的交点 (作网格时可要求界点为两族网格交点). 如图 2.6, $P_0(x_{i_0}, y_{j_0})$ 是界点, $P_1(x_{i_0 + 1}, y_{j_0})$ 和 $P_2(x_{i_0}, y_{j_0 - 1})$ 是与之相邻的内点. 过 $(x_{i_0 + \frac{1}{2}}, y_{j_0}), (x_{i_0}, y_{j_0 - \frac{1}{2}})$ 分别作与 $y$ 轴和 $x$ 轴平行的直线, 它们与外边界 $\Gamma$ 一起截出一曲边三角形 $\triangle ABC$ , 其为相应于 $P_0$ 的对偶单元. 于 $\triangle ABC$ 积分 (2.53) 两端, 并利用 Green 公式, 得

$$
- \int_ {\widehat {A B C A}} \frac {\partial u}{\partial \boldsymbol {n}} d s = \iint_ {\widetilde {\Delta} _ {A B C}} f d x d y. \tag {2.71}
$$

![](images/ca89028d5128624222dcbc38d413d4eebf0e8ca8b3115383e52f9a3f2c8630b4.jpg)  
图2.6 相应于 $P_0$ 的对偶单元

而

$$
\begin{array}{l} \int_ {\overline {{A B}}} \frac {\partial u}{\partial \boldsymbol {n}} d s \approx \frac {u _ {P _ {2}} - u _ {P _ {0}}}{h _ {2}} \cdot \overline {{A B}}, \\ \int_ {\overline {{B C}}} \frac {\partial u}{\partial \boldsymbol {n}} d s \approx \frac {u _ {P _ {1}} - u _ {P _ {0}}}{h _ {1}} \cdot \overline {{B C}}, \\ \int_ {\widehat {C A}} \frac {\partial u}{\partial \boldsymbol {n}} d s = \int_ {\widehat {C A}} (\gamma - k u) d s \approx \left(\gamma_ {P _ {0}} - k _ {P _ {0}} u _ {P _ {0}}\right) \cdot \widehat {C A}, \\ \end{array}
$$

以之代到 (2.71) 即得逼近 (2.70) 的差分方程:

$$
- \left[ \frac {u _ {p _ {2}} - u _ {p _ {0}}}{h _ {2}} \cdot \overline {{A B}} + \frac {u _ {p _ {1}} - u _ {p _ {0}}}{h _ {1}} \cdot \overline {{B C}} + \left(\gamma_ {p _ {0}} - k _ {p _ {0}} u _ {p _ {0}}\right) \cdot \widehat {C A} \right] = \iint_ {\widetilde {\Delta} A B C} f d x d y. \tag {2.72}
$$

可见用有限体积法处理第二、第三边值条件特别方便. 对非正则内点, 仍可像前面一样建立一形如 (2.68) 或 (2.69) 的方程

# 2.3.3 习题

1. 用有限体积法构造逼近方程

$$
- \nabla (k \nabla u) = - \left[ \frac {\partial}{\partial x} \left(k \frac {\partial u}{\partial x}\right) + \frac {\partial}{\partial y} \left(k \frac {\partial u}{\partial y}\right) \right] = f \tag {2.73}
$$

的第一边值问题的五点差分格式, 这里 $k = k(x,y) \geqslant k_{\min} > 0$

2. 用有限体积法构造逼近方程 (2.73) 的第二边值问题的五点差分格式  
3. 设步长 $h_1 = h_2 = h$ . 记

$$
\begin{array}{l} \diamond u _ {i j} = u _ {i + 1, j} + u _ {i, j + 1} + u _ {i - 1, j} + u _ {i, j - 1} - 4 u _ {i j}, \\ \square u _ {i j} = u _ {i + 1, j + 1} + u _ {i - 1, j + 1} + u _ {i - 1, j - 1} + u _ {i + 1, j - 1} - 4 u _ {i j}. \\ \end{array}
$$

证明逼近 Laplace 方程 $\Delta u = 0$ 的差分方程

$$
\frac {1}{6 h ^ {2}} \left(4 \diamond u _ {i, j} + \square u _ {i, j}\right) = 0
$$

的截断误差的绝对值

$$
\left| R _ {i, j} (u) \right| = \frac {4 0 h ^ {6}}{3 \cdot 8 !} \theta M _ {8},
$$

其中 $M_8$ 是 $u$ 的八阶偏导数的绝对值于考虑区域的上确界， $|\theta |\leqslant 1$

# 2.4 三角网的差分格式

利用有限体积法, 可将矩形网的差分格式推广到三角网, 得到三角网的差分格式, 文献上也称之为广义差分法 (见 [9, 20]). 三角网的差分格式具有网格灵活、边值条件容易处理等优点, 特别地, 它还保持积分守恒 (质量守恒), 所以受到应用部门的欢迎.

考虑有界域 $G$ 上的 Poisson 方程

$$
- \Delta u = f. \tag {2.74}
$$

在边界 $\Gamma$ 上满足第一、第二或第三边值条件. 我们先对 $G$ 作三角剖分. 如图2.7, 在 $\Gamma$ 上取一系列节点, 以它们为顶点作逼近 $\Gamma$ 的闭折线 $\widetilde{\Gamma}$ . 设 $\widetilde{G}$ 为由 $\widetilde{\Gamma}$ 围成的逼近 $G$ 的多边形域, 然后把 $\widetilde{G}$ 分割成有限个小三角形之和, 使不同三角形无重叠的内部区域, 任一三角形的顶点不属于其他三角形边的内部. 此外, 我们还要求这些三角形的内角不大于 $90^{\circ}$ . 这样, 就把 $\widetilde{G}$ 分割成一三角网, 称为 $G$ 的三角剖分. 组成剖分的小三角形也称为三角单元.

![](images/17d5eb1ed62a01297896f630bfb16c88b3b9c1a94a50da5068a5b3e50fa9f93b.jpg)  
图2.7 三角剖分

对于任一节点, 考虑所有以它为顶点的三角单元和以它为端点的三角形的边. 过每边作中垂线, 它们依次交于相应三角形的外心, 从而得到围绕该节点的小多边形域, 称为对偶单元 (也叫控制体积). 这些对偶单元全体构成区域 $G$ 的一个新的网格剖分, 称为对偶剖分. 如图2.8(a), 内点 $P_{0} \in G$ 的对偶单元是六边形域, $P_{0}$ 在六边形内部. 若 $P_{0}$ 是界点, 即 $P_{0} \in \Gamma$ , 则 $P_{0}$ 是其对偶单元的一个顶点 (参看图2.8(b)).

现就每一内点建立差分方程.如图2.8(a)，设 $P_0$ 是内点， $P_{1},P_{2},\dots ,P_{6}$ 是和 $P_0$ 相邻的节点， $q_{i}$ 为三角形 $P_{0}P_{i}P_{i + 1}$ $(P_{7} = P_{1}$ ，下同）的外心， $m_i$ 是线段 $\overline{P_0P_i}$ 的中点， $K_{P_0}^*$ 是由六边形 $q_{1}q_{2}\dots q_{6}$ 围成的对偶单元.在子域 $K_{P_0}$ 积分(2.74)，得

$$
- \iint_ {K _ {P _ {0}} ^ {*}} \Delta u \mathrm {d} x \mathrm {d} y = \iint_ {K _ {P _ {0}} ^ {*}} f \mathrm {d} x \mathrm {d} y.
$$

![](images/67bd657eadec09299939e7a383173dd800060ff42f81aef72d5076702616078d.jpg)  
(a) 内部对偶

![](images/296d09290ab75762beb3955dd3b8db95478acc72fc75c5eb017cf00d701201ad.jpg)  
(b) 边界对偶单元  
图2.8 对偶单元

利用Green公式，可将上式改写成

$$
- \int_ {\partial K _ {P _ {0}} ^ {*}} \frac {\partial u}{\partial \boldsymbol {n}} \mathrm {d} s = \iint_ {K _ {P _ {0}} ^ {*}} f \mathrm {d} x \mathrm {d} y. \tag {2.75}
$$

$\partial K_{P_0}^*$ 是 $K_{P_0}^*$ 的边界, $\pmb{n}$ 是 $\partial K_{P_0}$ 的单位外法向量. 注意

$$
\int_ {\partial K _ {P _ {0}} ^ {*}} \frac {\partial u}{\partial \boldsymbol {n}} d s = \sum_ {i = 1} ^ {6} \int_ {\overline {{q _ {i} q _ {i + 1}}}} \frac {\partial u}{\partial \boldsymbol {n}} d s = \sum_ {i = 1} ^ {6} \frac {\overline {{q _ {i} q _ {i + 1}}}}{P _ {0} P _ {i + 1}}. \tag {2.76}
$$

$$
\left[ u \left(P _ {i + 1}\right) - u \left(P _ {0}\right) \right] + m \left(G _ {0}\right) R _ {G _ {0}} (u) \quad \left(q _ {7} = q _ {1}\right).
$$

以之代到 (2.75)，即得点 $P_0$ 的差分方程：

$$
- \sum_ {i} \frac {\overline {{q _ {i} q _ {i + 1}}}}{\overline {{P _ {0} P _ {i + 1}}}} \left(u _ {P _ {i + 1}} - u _ {P _ {0}}\right) = \iint_ {K _ {P _ {0}} ^ {*}} f d x d y = m \left(K _ {P _ {0}} ^ {*}\right) \cdot \varphi_ {0}, \tag {2.77}
$$

$$
\varphi_ {0} = \frac {1}{m \left(K _ {P _ {0}} ^ {*}\right)} \iint_ {K _ {P _ {0}} ^ {*}} f d x d y,
$$

其中 $m(K_{P_0}^*)$ 是 $K_{P_0}^*$ 的面积.其次建立界点的差分方程.如图2.8(b)，设 $P_0$ 是界点，相应的对偶单元为由 $P_0m_1q_1q_2q_3m_4P_0$ 围成的多边形 $K_{P_0}^*$ .若在 $\Gamma$ 上给的是第一边值条件(2.54)，则令

$$
u _ {P _ {0}} = \alpha \left(P _ {0}\right). \tag {2.78}
$$

若给的是第二或第三边值条件，例如

$$
\left. \frac {\partial u}{\partial n} + k u \right| _ {r} = \gamma \tag {2.79}
$$

$(k \equiv 0$ 就是第二边值条件), 则需补充一个方程. 如图 2.8(b), 此时和 (2.75) 类似地有

$$
\begin{array}{l} - \left(\int_ {\overline {{m _ {4} P _ {0}}}} \frac {\partial u}{\partial \boldsymbol {n}} d s + \int_ {\overline {{P _ {0} m _ {1}}}} \frac {\partial u}{\partial \boldsymbol {n}} d s + \int_ {\overline {{m _ {1} q _ {1}}}} \frac {\partial u}{\partial \boldsymbol {n}} d s + \right. \\ \left. \int_ {\overline {{q _ {1} q _ {2}}}} \frac {\partial u}{\partial \boldsymbol {n}} \mathrm {d} s + \int_ {\overline {{q _ {2} q _ {3}}}} \frac {\partial u}{\partial \boldsymbol {n}} \mathrm {d} s + \int_ {\overline {{q _ {3} m _ {4}}}} \frac {\partial u}{\partial \boldsymbol {n}} \mathrm {d} s\right) \\ = \iint_ {K _ {P _ {0}} ^ {*}} f \mathrm {d} x \mathrm {d} y. \tag {2.80} \\ \end{array}
$$

上式左端括号内后四项仿照公式（2.76）的方法离散化，例如

$$
\int_ {\overline {{m _ {1} q _ {1}}}} \frac {\partial u}{\partial \boldsymbol {n}} d s \approx \frac {\overline {{m _ {1} q _ {1}}}}{\overline {{P _ {0} P _ {1}}}} (u _ {1} - u _ {0}), \quad \int_ {\overline {{q _ {1} q _ {2}}}} \frac {\partial u}{\partial \boldsymbol {n}} d s \approx \frac {\overline {{q _ {1} q _ {2}}}}{\overline {{P _ {0} P _ {2}}}} (u _ {2} - u _ {0}).
$$

(2.80) 左端前两项是沿外边界的积分, 利用条件 (2.79) 消去法向导数, 得到仅含 $u$ 的积分. 假定 $k, \gamma$ 是常数, $u$ 用三角单元边上的线性函数去逼近, 则可利用梯形公式计算相应的积分, 于是得

$$
\begin{array}{l} \int_ {\overline {{m _ {4} P _ {0}}}} \frac {\partial u}{\partial \boldsymbol {n}} d s = \int_ {\overline {{m _ {4} P _ {0}}}} (\gamma - \dot {k} u) d s \\ = \overline {{m _ {4} P _ {0}}} \left[ \gamma - \frac {1}{2} k \left(u _ {P _ {0}} + u _ {m _ {4}}\right) \right] \\ = \frac {1}{2} \overline {{P _ {0} P _ {4}}} \left[ \gamma - \frac {1}{4} k \left(3 u _ {P _ {0}} + u _ {P _ {4}}\right) \right], \tag {2.81} \\ \end{array}
$$

同理

$$
\int_ {\overline {{P _ {0} m _ {1}}}} \frac {\partial u}{\partial \boldsymbol {n}} d s = \frac {1}{2} \overline {{P _ {0} P _ {1}}} \left[ \gamma - \frac {1}{4} k \left(3 u _ {P _ {0}} + u _ {P _ {1}}\right) \right]. \tag {2.82}
$$

将这些公式代到 (2.80)，就得到界点的差分方程。显然所有内点、界点的差分方程组成一个封闭的线性代数方程组，其系数矩阵是对称的稀疏矩阵。

从前面的推导看出，用有限体积法构造三角网差分格式和矩形网情形完全类似，而边值条件的处理则更方便灵活，特别是第二、第三边值条件，其处理方法跟内点没有实质区别。还应指出，方程（2.75）是（2.74）的积分形式，表示某一物理量在单元 $G_{0}$ 守恒，而差分方程（2.77）则是守恒律的离散形式。

# 例2.1 五点差分格式

在矩形网上用同向对角线将每一矩形单元分成两个直角三角形，得到“直角三角剖分”，对偶单元是矩形，由此导出的差分格式恰是五点差分格式.

# 例2.2 正三角网上的差分格式

图2.9是一个正三角网，每个三角形的边长为 $h$ 。取一内节点 $P_0$ ，设与之相邻的六个

节点为 $P_{1}, P_{2}, P_{3}, P_{4}, P_{5}, P_{6}$ . 过 $\overline{P_{0}P_{i}}$ 的中点 $m_{i}$ 作 $\overline{P_{0}P_{i}}$ 的中垂线, 依次交于 $\triangle P_{0}P_{i}P_{i+1}$ 的外心 $q_{i}$ (这时与重心重合). 正六边形 $q_{1}q_{2}\cdots q_{6}$ 围成的多边形域 $G_{0}$ 是围绕 $P_{0}$ 的对偶单元. 显然

$$
\overline {{P _ {0} P _ {i}}} = h, \quad \overline {{q _ {i} q _ {i + 1}}} = \frac {h}{\sqrt {3}}.
$$

$G_{0}$ 的面积 $m\left(G_{0}\right) = \frac{\sqrt{3}h^{2}}{2}$ . 因此差分格式 (2.77) 为

$$
- \frac {2}{3 h ^ {2}} \left(\sum_ {i = 1} ^ {6} u _ {i} - 6 u _ {0}\right) = \frac {2}{\sqrt {3} h ^ {2}} \iint_ {G _ {0}} f (x, y) \mathrm {d} x \mathrm {d} y. \tag {2.83}
$$

![](images/1e61297fc3e75349338fc611d2fcfb5d544b4c15a9db5e40ac43cb6cb4164ff9.jpg)  
图2.9 正三角网

例2.3 正六边形网上的差分格式

图2.10是正六边形网，边长为 $h$ 设 $P_0$ 是任一内节点，与之相邻的三个节点为 $P_{1},P_{2},P_{3}$ 过 $\overline{P_0P_i}$ 的中点作中垂线，彼此相交于六边形的中心 $q_{1},q_{2},q_{3},\triangle q_{1}q_{2}q_{3}$ 是围绕 $P_0$ 的对偶单元.显然 $\overline{P_0P_i} = h,\overline{q_iq_{i + 1}} = \sqrt{3} h,$ 对偶单元 $\triangle q_{1}q_{2}q_{3}$ 的面积为 $\frac{3\sqrt{3}h^2}{4}$

![](images/77af198dbb59c9dd1e43b5fa222289f9055901d8f5e38a6ce4b2f085e1295b96.jpg)  
图2.10 正六边形网

因此差分格式为

$$
- \frac {4}{3 h ^ {2}} \left(\sum_ {i = 1} ^ {3} u _ {i} - 3 u _ {0}\right) = \frac {4}{3 \sqrt {3} h ^ {2}} \iint_ {\triangle q _ {1} q _ {2} q _ {3}} f (x, y) \mathrm {d} x \mathrm {d} y. \tag {2.84}
$$

# 2.4.1 习题

1. 试证三角网的差分格式（第一或第三边值条件）的系数矩阵对称  
2.构造逼近（2.73）的三角网的差分格式  
3. 求出差分格式 (2.83) 和 (2.84) 的截断误差的阶 (分别为 $O(h^{2})$ 和 $O(h)$ ).

# *2.5 极值定理和敛速估计

为了得到差分解的收敛性、敛速估计及其稳定性，需要对差分解作某种先验估计，极值定理是作这类估计的常用方法.

# 2.5.1 差分方程

考虑二阶椭圆偏微分方程第一边值问题：

$$
\left\{ \begin{array}{l} - \left(A u _ {x} ^ {\prime}\right) _ {x} - \left(B u _ {y} ^ {\prime}\right) _ {y} + C u _ {x} ^ {\prime} + D u _ {y} ^ {\prime} + E u = F, \quad (x, y) \in G, \\ u | _ {\Gamma} = \alpha , \end{array} \right. \tag {2.85}
$$

其中 $A(x,y),B(x,y)$ 属于 $C^1 (\bar{G})$ ， $C(x,y),D(x,y),E(x,y)$ 和 $F(x,y)$ 属于 $C(\bar{G})$ ， $\alpha \in C(\Gamma)$ ，且 $A(x,y)\geqslant A_{\min} > 0,B(x,y)\geqslant B_{\min} > 0,E\geqslant 0.$

如2.3节构造矩形网， $h_1$ 和 $h_2$ 分别为沿 $x$ 和 $y$ 方向的步长.用 $G_{h}$ 表示网格内点集合， $\Gamma_h$ 表示网格界点集合， $\bar{G}_h = G_h\cup \Gamma_h$ .本节总假定 $G_{h}$ 是连通的，就是说，对任意两节点 $\bar{P},\bar{\bar{P}}\in G_{h}$ ，必有一串节点 $P_{i}\in G_{h}(i = 1,2,\dots ,m - 1)$ ，可与 $\bar{P},\bar{\bar{P}}$ 排成下列顺序：

$$
\bar {P}, P _ {1}, P _ {2}, \dots , P _ {m - 1}, \bar {\bar {P}},
$$

使前后两点为相邻节点

对于正则内点 $(x_{i},y_{j})$ ，用如下的差分方程逼近(2.85)：

$$
- \left[ A _ {i - \frac {1}{2}, j} (u _ {i j}) _ {\bar {x}} \right] _ {x} - \left[ B _ {i, j - \frac {1}{2}} (u _ {i j}) _ {\bar {y}} \right] _ {y} + C _ {i j} (u _ {i j}) _ {\widehat {x}} + D _ {i j} (u _ {i j}) _ {\widehat {y}} + E _ {i j} u _ {i j} = F _ {i j}, \tag {2.86}
$$

其中

$$
A _ {i - \frac {1}{2}, j} = A \left(x _ {i - \frac {1}{2}}, y _ {j}\right) = A \left(\left(i - \frac {1}{2}\right) h _ {1}, j h _ {2}\right),
$$

$$
B _ {i, j - \frac {1}{2}} = B \left(x _ {i}, y _ {j - \frac {1}{2}}\right) = B \left(i h _ {1}, \left(j - \frac {1}{2}\right) h _ {2}\right),
$$

而

$$
\left\{ \begin{array}{l} \left(u _ {i j}\right) _ {\bar {x}} = \frac {u _ {i j} - u _ {i - 1 , j}}{h _ {1}}, \quad \left(u _ {i j}\right) _ {\bar {y}} = \frac {u _ {i j} - u _ {i , j - 1}}{h _ {2}}, \\ \left(u _ {i j}\right) _ {x} = \frac {u _ {i + 1 , j} - u _ {i , j}}{h _ {1}}, \quad \left(u _ {i j}\right) _ {y} = \frac {u _ {i , j + 1} - u _ {i , j}}{h _ {2}}, \\ \left(u _ {i j}\right) _ {\bar {x}} = \frac {u _ {i + 1 , j} - u _ {i - 1 , j}}{2 h _ {1}}, \quad \left(u _ {i j}\right) _ {\bar {y}} = \frac {u _ {i , j + 1} - u _ {i , j - 1}}{2 h _ {2}}. \end{array} \right. \tag {2.87}
$$

显然，截断误差的阶为 $O\left(h_1^2 + h_2^2\right)$ 。方程（2.86）可改写为

$$
a _ {i j} u _ {i j} - \left(a _ {i - 1, j} u _ {i - 1, j} + a _ {i, j - 1} u _ {i, j - 1} + a _ {i + 1, j} u _ {i + 1, j} + a _ {i, j + 1} u _ {i, j + 1}\right) = F _ {i j}, \tag {2.88}
$$

其中

$$
\left\{ \begin{array}{l} a _ {i - 1, j} = h _ {1} ^ {- 2} \left(A _ {i - \frac {1}{2}, j} + \frac {h _ {1}}{2} C _ {i j}\right), \\ a _ {i, j - 1} = h _ {2} ^ {- 2} \left(B _ {i, j - \frac {1}{2}} + \frac {h _ {2}}{2} D _ {i j}\right), \\ a _ {i + 1, j} = h _ {1} ^ {- 2} \left(A _ {i + \frac {1}{2}, j} - \frac {h _ {1}}{2} C _ {i j}\right), \\ a _ {i, j + 1} = h _ {2} ^ {- 2} \left(B _ {i, j + \frac {1}{2}} - \frac {h _ {2}}{2} D _ {i j}\right), \\ a _ {i j} = h _ {1} ^ {- 2} \left(A _ {i + \frac {1}{2}, j} + A _ {i - \frac {1}{2}, j}\right) + h _ {2} ^ {- 2} \left(B _ {i, j + \frac {1}{2}} + B _ {i, j - \frac {1}{2}}\right) + E _ {i j}. \end{array} \right. \tag {2.89}
$$

由系数 $A, B$ 的假设条件，只要 $h_1$ 和 $h_2$ 充分小，则 $a_{i-1,j}, a_{i,j-1}, a_{i+1,j}, a_{i,j+1}$ 和 $a_{ij}$ 均大于0，且

$$
a _ {i j} - \left(a _ {i - 1, j} + a _ {i, j - 1} + a _ {i + 1, j} + a _ {i, j + 1}\right) = E _ {i j} \geqslant 0. \tag {2.90}
$$

对于非正则内点, 则建立一不等距差分方程. 例如设 $(x_{i},y_{j})$ 为图2.5中点“0”, 则用

$$
\left[ A _ {i - \frac {1}{2}, j} (u _ {i j}) _ {\bar {x}} \right] _ {x} = \frac {1}{\bar {h} _ {1}} \left(A _ {i + \frac {1}{2}, j} \frac {u _ {i + 1 , j} - u _ {i j}}{h _ {1}} - A _ {i - \frac {1}{2}, j} \frac {u _ {i , j} - u _ {i - 1 , j}}{h _ {1} ^ {-}}\right)
$$

和

$$
(u _ {i j}) _ {\widehat {x}} = \frac {u _ {i + 1 , j} - u _ {i - 1 , j}}{h _ {1} + h _ {1} ^ {-}}
$$

依次代替 (2.86) 中的相应项, 其中 $h_1^-, \bar{h}_1$ 如 (2.68). 此时仍可将 (2.86) 写成形式 (2.88). 只要 $h_1, h_1^-, h_2, h_2^-$ 充分小, 则 (2.88) 左端之系数 $a_{i-1,j}, a_{i,j-1}, a_{i+1,j}, a_{i,j+1}$ 和 $a_{ij}$ 就是正的, 且 (2.90) 成立. 显然, 在非正则内点, 差分逼近的阶为 $O(h_1 + h_2)$ .

在研究 (2.88) 的极值性质之前, 我们指出它的几个简单而有用的性质. 将 (2.88) 改写成

$$
L _ {h} u _ {i j} = a _ {i j} u _ {i j} - a _ {i - 1, j} u _ {i - 1, j} - a _ {i, j - 1} u _ {i, j - 1} - a _ {i + 1, j} u _ {i + 1, j} - a _ {i, j + 1} u _ {i, j + 1} = F _ {i j}. \tag {2.91}
$$

将网格内点按适当次序排列, 例如从左下角网点开始, 按由左向右、由下向上的顺序排列, 得到一线性代数方程组, 其系数矩阵 $\mathbf{A}$ 有下列性质:

(i) $A$ 的每行最多有五个非零元素, 所以 $A$ 为稀疏矩阵. 在排列网点顺序时, 应尽量使非零元素“靠近”对角线, 这对消元法特别有利.  
(ii) $\mathbf{A}$ 的对角元素是正的, 非对角元素是非正的. 非对角元素绝对值之和不超过对角元素, 即

$$
a _ {i - 1, j} + a _ {i, j - 1} + a _ {i + 1, j} + a _ {i, j + 1} \leqslant a _ {i j}. \tag {2.92}
$$

当 $(i,j)$ 为非正则内点时，其四个相邻点至少有一个是界点，比如设 $(i - 1,j)$ 是界点，则可将(2.91）中之相应项 $a_{i - 1,j}u_{i - 1,j}$ 移到右端，视(2.91）左端之 $a_{i - 1,j} = 0$ ，此时不等式(2.92）严格小于号成立.所以矩阵 $\mathbf{A}$ 对角占优，这对于保证迭代法的收敛性是重要的.

(iii) 若方程 (2.85) 对称, 即 $C = D = 0$ , 则矩阵 $\mathbf{A}$ 也对称（设非正则内点的格式为 (2.69)).

上述性质在构造差分方程的解法时特别有用

# 2.5.2 极值定理

现在讨论差分方程 (2.91) 的极值性质

定理2.2（极值定理）设 $u_{ij}$ 是 $\bar{G}_h$ 上的任一网格函数.若 $L_{h}u_{ij}\leqslant 0(L_{h}u_{ij}\geqslant 0)$ 对任意 $(x_i,y_j)\in G_h$ ，则 $u_{ij}$ 不可能在内点取正的极大值(负的极小值)，除非 $u_{ij}\equiv$ 常数

证明 只证明定理的第一部分, 因为第二部分是类似的. 用反证法, 设 $u_{ij} \neq$ 常数, $u_{ij}$ 在 $G_h$ 中某点达到正的极大值 $M$ . 由于 $G_h$ 连通, 必有某一内点 $(x_{i_0}, y_{j_0})$ , 使 $u_{i_0j_0} = M$ , 且至少有一个相邻网点, 比如 $(x_{i_0 - 1}, y_{j_0})$ , 使 $u_{i_0 - 1,j_0} < M$ . 于是

$$
L _ {h} u _ {i _ {0} j _ {0}} > \left(a _ {i _ {0} j _ {0}} - a _ {i _ {0} - 1, j _ {0}} - a _ {i _ {0} + 1, j _ {0} - 1} - a _ {i _ {0} + 1, j _ {0}} - a _ {i _ {0}, j _ {0} + 1}\right) M \geqslant 0
$$

(若 $a_{i_0 - 1,j_0} = 0$ ，则 $\geqslant$ 改为 $\geqslant ,\geqslant$ 改为 $\succ$ )与假设矛盾

推论2.1 差分方程(2.91）有唯一解

证明只需证明相应的齐问题(边值和右端都恒等于0)只有平凡解.实际上，设 $u_{ij}$ 是齐问题的解，则由定理2.2, $u_{ij}$ 既不能在 $G_{h}$ 取正的极大值，也不能取负的极小值，因此 $u_{ij} \equiv 0$

推论2.2 若网格函数 $u_{ij}$ 满足

$$
L _ {h} u _ {i j} \geqslant 0, \quad \forall (x _ {i}, y _ {j}) \in G _ {h},
$$

$$
u _ {i j} \geqslant 0, \quad \forall (x _ {i}, y _ {j}) \in \Gamma_ {h},
$$

则 $u_{ij} \geqslant 0, \forall (x_i, y_j) \in G_h$ .

证明 由定理2.2直接推得

定理2.3（比较定理）设 $u_{ij}$ 和 $U_{ij}$ 是两个网格函数，满足

$$
\left| L _ {h} u _ {i j} \right| \leqslant L _ {h} U _ {i j}, \quad \forall (x _ {i}, y _ {j}) \in G _ {h}, \tag {2.93}
$$

$$
\left| u _ {i j} \right| \leqslant U _ {i j}, \quad \forall (x _ {i}, y _ {j}) \in \Gamma_ {h}, \tag {2.94}
$$

则

$$
\left| u _ {i j} \right| \leqslant U _ {i j}, \quad \forall \left(x _ {i}, y _ {j}\right) \in G _ {h}. \tag {2.95}
$$

证明 由(2.93)和(2.94)可知

$$
\left\{ \begin{array}{l l} {L _ {h} \left(U _ {i j} - u _ {i j}\right) \geqslant 0,} & {\text {于} G _ {h},} \\ {U _ {i j} - u _ {i j} \geqslant 0,} & {\text {于} \varGamma_ {h}} \end{array} \right.
$$

和

$$
\left\{ \begin{array}{l l} {L _ {h} \left(U _ {i j} + u _ {i j}\right) \geqslant 0,} & {\text {于} G _ {h},} \\ {U _ {i j} + u _ {i j} \geqslant 0,} & {\text {于} \Gamma_ {h}.} \end{array} \right.
$$

由推论2.2便知(2.95)成立

推论2.3 差分方程

$$
\left\{ \begin{array}{l l} {L _ {h} u _ {i j} = 0,} & {\text {于} G _ {h},} \\ {u _ {i j} = \alpha_ {i j},} & {\text {于} \Gamma_ {h}} \end{array} \right.
$$

的解 $u_{ij}$ 满足不等式

$$
\max  _ {G _ {h}} | u _ {i j} | \leqslant \max  _ {\Gamma_ {h}} | \alpha_ {i j} |. \tag {2.96}
$$

证明 设 $U_{ij}$ 是下列问题的解：

$$
L _ {h} U _ {i j} = 0, \quad \text {于} G _ {h},
$$

$$
U _ {i j} = | \alpha_ {i j} |  , \quad \text {于}   \Gamma_ {h}  ,
$$

则由定理2.3

$$
| u _ {i j} | \leqslant U _ {i j},   \text {于}   \bar {G} _ {h}.
$$

若 $U_{ij} \equiv$ 常数, 则 $U_{ij} \equiv \max_{\Gamma_h} |\alpha_{ij}|$ . 若 $U_{ij} \neq$ 常数, 则由定理 2.2, 函数 $U_{ij} (\geqslant 0)$ 的最大值只能在 $\Gamma_h$ 达到, 因此 $U_{ij} \leqslant \max_{\Gamma_h} |\alpha_{ij}|$ , 从而 (2.96) 成立.

# 2.5.3 五点差分格式的敛速估计

设 $u = u(x,y)$ 是Poisson方程第一边值问题的解（参看(2.53)(2.54)), $u_{ij}$ 是五点差分格式

$$
\left\{ \begin{array}{l l} {L _ {h} u _ {i j} = \varphi_ {i j},} & {\text {于} G _ {h},} \\ {u _ {i j} = \alpha_ {i j},} & {\text {于} \Gamma_ {h}} \end{array} \right.
$$

的解, 此时 $A = B = 1$ , $C = D = E = 0$ . 当 $(x_i, y_j)$ 是正则内点, 且 $h_1 = h_2$ 时, (2.91) 中的 $a_{ij} = 4$ , $a_{i-1,j} = a_{i,j-1} = a_{i+1,j} = a_{i,j+1} = 1$ . 当 $(x_i, y_j)$ 为非正则内点时, 需对系数作适当修正. 令 $e_{ij} = u(x_i, y_j) - u_{ij}$ , 若 $u \in C^4(\bar{G})$ , 则 $e_{ij}$ 满足

$$
\left\{ \begin{array}{l l} {L _ {h} e _ {i j} = R _ {i j},} & {\text {于} G _ {h},} \\ {e _ {i j} = 0,} & {\text {于} \Gamma_ {h}.} \end{array} \right. \tag {2.97}
$$

截断误差

$$
R _ {i j} = \left\{ \begin{array}{l l} {O \left(h _ {1} ^ {2} + h _ {2} ^ {2}\right),} & {(x _ {i}, y _ {j}) \text {是 正 则 内 点},} \\ {O \left(h _ {1} + h _ {2}\right),} & {(x _ {i}, y _ {j}) \text {是 非 正 则 内 点}.} \end{array} \right.
$$

令 $h = \sqrt{h_1^2 + h_2^2}$ ，可设

$$
\left| R _ {i j} \right| \leqslant K h,
$$

其中 $K$ 是常数

现在估计误差 $e_{ij}$ .不妨设 $(0,0)\in G,R$ 是以 $(0,0)$ 为圆心且包含 $G$ 在内部的最小圆域的半径.令

$$
E _ {i j} = \frac {K h}{4} \left(R ^ {2} - x _ {i} ^ {2} - y _ {j} ^ {2}\right).
$$

显然, $E_{ij}$ 于 $\overline{G_h}$ 非负. 又 $x_i^2$ 关于 $x$ 方向的二阶中心差商等于2, $y_j^2$ 关于 $y$ 方向的二阶中心差商也等于2, 故当 $(x_i,y_j)$ 是正则内点时, $L_{h}E_{ij} = Kh$ . 而当 $(x_i,y_j)$ 为非正则内点时, $L_{h}$ 中仍出现二阶差商 (但不一定是中心差商), 此时仍有 $L_{h}E_{ij} = Kh$ . 这样, $E_{ij}$ 满足

$$
\left\{ \begin{array}{l l} L _ {h} E _ {i j} = K h, & \left(x _ {i}, y _ {j}\right) \in G _ {h}, \\ E _ {i j} \geqslant 0, & \left(x _ {i}, y _ {j}\right) \in \Gamma_ {h}. \end{array} \right. \tag {2.98}
$$

将定理2.3用于(2.97)和(2.98)，便得

$$
\left| e _ {i j} \right| \leqslant E _ {i j} = \frac {K h}{4} \left(R ^ {2} - x _ {i} ^ {2} - y _ {j} ^ {2}\right),
$$

从而

$$
\max  _ {G _ {h}} \left| e _ {i j} \right| \leqslant \frac {K R ^ {2} h}{4}. \tag {2.99}
$$

足见若 $u(x,y)\in C^4 (\bar{G})$ ，则差分解 $u_{ij}$ 一致收敛到 $\pmb{u}$ 且有敛速估计(2.99).

注2.3 极值定理在证明椭圆型差分方程的稳定性及差分解的误差估计中起重要作用. 在下列习题1,2中给出极值定理的其他形式.

# 2.5.4 习题

1. 设 $\bar{I}_h = \{x_i : i = 0, 1, \dots, N, x_0 < x_1 < \dots < x_N\}$ , $y_i$ 是 $\bar{I}_h$ 上的网格函数. 又

$$
l y _ {i} = - \left(a _ {i} y _ {i - 1} - b _ {i} y _ {i} + c _ {i} y _ {i + 1}\right) + q _ {i} y _ {i}, i = 1, 2, \dots , N - 1,
$$

其中 $a_{i}, b_{i}, c_{i}$ 恒正, $q_{i}$ 非负, 且 $a_{i} + c_{i} \leqslant b_{i}$ . 证明当 $ly_{i} \leqslant 0 (ly_{i} \geqslant 0)$ 时, $y_{i}$ 不能在内点取正的极大值 (负的极小值), 除非 $y_{i}$ 等于常数.

2. 在上题中, 若设 $d_{i} = b_{i} - a_{i} - c_{i} + q_{i} > 0 (i = 1,2,\dots ,N - 1)$ , 则差分方程

$$
\left\{ \begin{array}{l l} l y _ {i} = \varphi_ {i}, & i = 1, 2, \dots , N - 1, \\ y _ {0} = y _ {N} = 0 & \end{array} \right.
$$

的解满足

$$
\max  _ {i} | y _ {i} | \leqslant \max  _ {i} \frac {| \varphi_ {i} |}{d _ {i}}.
$$

3. 利用上题估计差分方程 (2.28) (2.29) 解的收敛阶 (假定 $r = 0, q \geqslant q_0 > 0, h_i \equiv h$ ).

# 第3章

# 抛物型方程的有限差分法

椭圆型方程描写的是状态（如温度、电位等）不随时间 $t$ 改变的问题，称为驻定问题。现在讨论与时间 $t$ 有关的非驻定问题。驻定问题可看成是某一非驻定问题当 $t \to \infty$ 时的渐近状态，所以当我们用渐近方法（例如迭代法）求解驻定问题时，只关心最终状态，而不管中间过程。相反，非驻定问题的瞬时状态有物理意义，需要我们求解。在考虑偏微分方程的数值解法时，注意到这两类问题的联系和区别是有益的。下面分别讨论抛物型方程和双曲型方程的差分法。

# 3.1 最简差分格式

考虑一维热传导方程

$$
\frac {\partial u}{\partial t} = a \frac {\partial^ {2} u}{\partial x ^ {2}} + f (x), \quad 0 <   t \leqslant T, \tag {3.1}
$$

其中 $a$ 是正常数， $f(x)$ 是给定的连续函数. 按照初边值条件的不同给法，可将 (3.1) 的定解问题分为两类：

第一，初值问题（也称Cauchy问题）：求具有一定阶偏微商的函数 $u(x,t)$ ，满足方程(3.1)和初值条件

$$
u (x, 0) = \varphi (x), \quad - \infty <   x <   \infty . \tag {3.2}
$$

第二，初边值问题（也称混合问题）：求具有一定阶偏微商的函数 $u(x,t)$ ，满足方程(3.1）和初值条件

$$
u (x, 0) = \varphi (x), \quad 0 <   x <   l \tag {3.3}
$$

及边值条件

$$
u (0, t) = u (l, t) = 0, \quad 0 \leqslant t \leqslant T. \tag {3.4}
$$

假定 $f(x)$ 和 $\varphi (x)$ 在相应区域光滑，并且在 $x = 0,x = l$ 相容 $(\varphi (0) = \varphi (l) = 0)$ ，则上述问题有唯一的光滑解.

现在考虑边值问题 (3.1) (3.3) (3.4) 的差分逼近. 取空间步长 $h = \frac{l}{J}$ 和时间步长 $\tau = \frac{T}{N}$ , 其中 $J, N$ 都是自然数. 用两族平行直线 $x = x_{j} = jh$ ( $j = 0,1,\dots,J$ ) 和 $t = t_{n} = n\tau (n = 0,1,\dots,N)$ 将矩形域 $\bar{G} = \{0 \leqslant x \leqslant l; 0 \leqslant t \leqslant T\}$ 分割成矩形网格, 网格节点为 $(x_{j}, t_{n})$ . 以 $G_{h}$ 表示网格内点集合, 即位于开矩形 $G$ 的网点集合; $\bar{G}_{h}$ 表示所有位于闭矩形 $\bar{G}$ 的网点集合; $\Gamma_{h} = \bar{G}_{h} - G_{h}$ 是网格界点集合 (参看图3.1).

其次，用 $u_{j}^{n}$ 表示定义在网点 $(x_{j},t_{n})$ 上的函数， $0\leqslant j\leqslant J,0\leqslant n\leqslant N.$ 用适当的差商代替方程（3.1）中相应的偏微商，便得到以下几种最简差分格式

![](images/1bc2291f20d9bfc777930b3ca328f254a1f8ba52989c72b0030048842e572720.jpg)  
图3.1 矩形网格

# （一）向前差分格式，即

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} + f _ {j}, \quad f _ {j} = f (x _ {j}), \tag {3.5}
$$

$$
u _ {j} ^ {0} = \varphi_ {j} = \varphi (x _ {j}), \quad u _ {0} ^ {n} = u _ {J} ^ {n} = 0, \tag {3.6}
$$

其中 $j = 1,2,\dots ,J - 1,n = 0,1,\dots ,N - 1.$ 以 $r = \frac{a\tau}{h^2}$ 表示网比.将(3.5)改写成便于计算的形式，使第 $n$ 层值（上标为 $n$ )在等式右边，第 $n + 1$ 层值在等式左边，则得

$$
u _ {j} ^ {n + 1} = r u _ {j + 1} ^ {n} + (1 - 2 r) u _ {j} ^ {n} + r u _ {j - 1} ^ {n} + \tau f _ {j}. \tag {3.7}
$$

取 $n = 0$ ，利用初值 $u_{j}^{0} = \varphi_{j}$ 和边值 $u_0^n = u_j^n = 0$ 可由(3.7)算出第一层值 $u_{j}^{1}$ .于(3.7)取 $n = 1$ ，又可利用 $u_{j}^{1}$ 和边值，由(3.7)算出 $u_{j}^{2}$ .如此下去，即可逐层求出所有 $u_{j}^{n}$ ，并视 $u_{j}^{n}$ 为精确解 $u(x_{j},t_{n})$ 的近似.由于第 $n + 1$ 层值通过第 $n$ 层值明显表为(3.7)，无需解线性代数方程组，如此的差分格式称为显格式.将(3.7)看成网点 $(x_{j},t_{n})$ 处的差分方程，它联系第 $n + 1$ 层的点 $(x_{j},t_{n + 1})$ 和第 $n$ 层的点 $(x_{j - 1},t_n),(x_j,t_n)$ 及 $(x_{j + 1},t_n)$ ，其分布如图3.2(a）所示.

记

$$
L u = \frac {\partial u}{\partial t} - a \frac {\partial^ {2} u}{\partial x ^ {2}},
$$

$$
L _ {h} ^ {(1)} u _ {j} ^ {n} = \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}.
$$

显然截断误差

$$
\begin{array}{l} R _ {j} ^ {n} (u) = L _ {h} ^ {(1)} u \left(x _ {j}, t _ {n}\right) - (L u) _ {j} ^ {n} \\ = - \tau \left(\frac {1}{1 2 r} - \frac {1}{2}\right) \left(\frac {\partial^ {2} \tilde {u}}{\partial t ^ {2}}\right) _ {j} ^ {n} + O \left(\tau^ {2} + h ^ {2}\right) \\ = O \left(\tau + h ^ {2}\right), \tag {3.8} \\ \end{array}
$$

其中 $\left(\frac{\partial^2\tilde{u}}{\partial t^2}\right)_j^n$ 是 $\frac{\partial^2u}{\partial t^2}$ 在矩形 $x_{j - 1} <   x <   x_{j + 1},t_n <   t <   t_{n + 1}$ 中的某点值

![](images/37414f366fbb6193c6b72fe12802a61608044b8ec0a60cfed7a2113017459174.jpg)  
(a) 向前差分格式

![](images/44e3d5f3f341550ff940e578bb31351cd21cda5f82ea1534caf11781177ad6d6.jpg)  
(b) 向后差分格式

![](images/4a2206503d4866187ea0735af442f4d24164de0f9007d650cecca01e44d422ce.jpg)  
(c)六点对称格式

![](images/eef0e8ca99d3ca4fb6c75c98a81d83f80a181df840f86968cf96112d3385e9e6.jpg)  
(d)Richardson格式   
图3.2 网点分布

# （二）向后差分格式，即

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} + f _ {j}, \tag {3.9a}
$$

$$
u _ {j} ^ {0} = \varphi_ {j} = \varphi (x _ {j}), \quad u _ {0} ^ {n} = u _ {J} ^ {n} = 0, \tag {3.9b}
$$

![](images/ddfa245586e0559666fc8792f880d4477c9a82f3a2e294f5d5bb88ab2d1f725e.jpg)

人物简介

其中 $j = 1,2,\dots ,J - 1,n = 0,1,\dots ,N - 1.$ 将(3.9a)改写为

$$
- r u _ {j + 1} ^ {n + 1} + (1 + 2 r) u _ {j} ^ {n + 1} - r u _ {j - 1} ^ {n + 1} = u _ {j} ^ {n} + \tau f _ {j}. \tag {3.10}
$$

令 $n = 0,1,2,\dots$ ，则可利用 $u_{j}^{0}$ 和边值确定 $u_{j}^{1}$ ，利用 $u_{j}^{1}$ 和边值确定 $u_{j}^{2}$ ，等等.现在第 $n + 1$ 层的值不能用第 $_n$ 层值明显表示，而是由线性代数方程组（3.10）确定，如此的差分格式称为隐格式.我们指出，(3.10)左端系数矩阵严格对角占优，方程总是可解的．将(3.9a）看成是网点 $(x_{j},t_{n + 1})$ 处的差分方程，它所联系的网点分布如图3.2(b）所示

令

$$
L _ {h} ^ {(2)} u _ {j} ^ {n} = \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}},
$$

则截断误差

$$
\begin{array}{l} R _ {j} ^ {n} (u) = L _ {h} ^ {(2)} u \left(x _ {j}, t _ {n}\right) - \left[ L u \right] _ {j} ^ {n} \\ = - \tau \left(\frac {1}{1 2 r} + \frac {1}{2}\right) \left(\frac {\partial^ {2} \tilde {u}}{\partial t ^ {2}}\right) _ {j} ^ {n} + O \left(\tau^ {2} + h ^ {2}\right) \\ = O \left(\tau + h ^ {2}\right), \tag {3.11} \\ \end{array}
$$

其中 $\left(\frac{\partial^2\tilde{u}}{\partial t^2}\right)_j^n$ 是 $\frac{\partial^2u}{\partial t^2}$ 在矩形 $x_{j - 1} <   x <   x_{j + 1},t_n <   t <   t_{n + 1}$ 中的某点值

（三）六点对称格式（Crank-Nicolson格式).将向前差分格式和向后差分格式作算术平均，即得六点对称格式：

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = \frac {a}{2} \left(\frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} + \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}\right) + f _ {j}, \tag {3.12a}
$$

$$
u _ {j} ^ {0} = \varphi_ {j} = \varphi (x _ {j}), \quad u _ {0} ^ {n} = u _ {J} ^ {n} = 0. \tag {3.12b}
$$

将 (3.12a) 改写为

$$
- \frac {r}{2} u _ {j + 1} ^ {n + 1} + (1 + r) u _ {j} ^ {n + 1} - \frac {r}{2} u _ {j - 1} ^ {n + 1} = \frac {r}{2} u _ {j + 1} ^ {n} + (1 - r) u _ {j} ^ {n} + \frac {r}{2} u _ {j - 1} ^ {n} + \tau f _ {j}, \tag {3.13}
$$

利用 $u_{j}^{0}$ 和边值便可逐层求到 $u_{j}^{n}$ . 六点对称格式是隐格式, 由第 $n$ 层计算第 $n + 1$ 层时, 需解线性代数方程组 (因系数矩阵严格对角占优, 方程组可唯一求解), 它所联系的网点分布如图3.2(c) 所示.

令

$$
L _ {h} ^ {(3)} u _ {j} ^ {n} = \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - \frac {a}{2} \left(\frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} + \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}\right).
$$

将截断误差

$$
R _ {j} ^ {n} (u) = L _ {h} ^ {(3)} u \left(x _ {j}, t _ {n}\right) - [ L u ] _ {j} ^ {n}
$$

于 $\left(x_{j},t_{n + \frac{1}{2}}\right)\left(t_{n + \frac{1}{2}} = \left(n + \frac{1}{2}\right)\tau\right)$ 展开,则得

$$
R _ {j} ^ {n} (u) = O \left(\tau^ {2} + h ^ {2}\right). \tag {3.14}
$$

（四）Richardson格式，即

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n - 1}}{2 \tau} = a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} + f _ {j}, \tag {3.15}
$$

或

$$
u _ {j} ^ {n + 1} = 2 r \left(u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}\right) + u _ {j} ^ {n - 1} + 2 \tau f _ {j}. \tag {3.16}
$$

这是三层显式差分格式，它联系的网点分布如图3.2(d）所示．显然截断误差的阶为 $O\left(\tau^{2} + h^{2}\right)$ .为了计算能够逐层进行，除初值 $u_{j}^{0}$ 外，还要用到 $u_{j}^{1}$ ，这可以用前述二层差分格式计算（为保证精度，可将 $[0,\tau ]$ 分成若干等份).

除以上四种差分格式外，还可作出（例如用待定系数法）许多逼近(3.1)，(3.3)一(3.4)解的差分格式，但并不是每一种差分格式都是可用的。衡量一个差分格式是否经济实用，由多方面的因素决定，主要有：

(1) 计算简单. 显格式无须解方程组, 计算较隐格式简单. 但隐格式 (3.10) 和 (3.13) 左端系数是与 $n$ 无关的三对角矩阵, 这相当于求解许多具不同右端但有同一三对角系数矩阵的方程组, 用消元法有许多方便. 再考虑到其他因素 (如稳定性), 隐格式也是可用的.  
(2) 收敛性和收敛速度. 当网比 $r$ 固定, 步长 $h \to 0$ 时, 差分解 $u_{j}^{n}$ 应收敛到精确解 $u(x_{j}, t_{n})$ , 并希望有尽可能快的收敛速度. 差分算子 $L_{h}$ 的截断误差的无穷小阶反映了

$L_{h}$ 对微分算子 $L$ 的逼近程度，因此可以期望，截断误差的阶越高，差分解的精度也越高.六点对称格式和Richardson格式截断误差的阶是 $O\left(\tau^{2} + h^{2}\right)$ ，而显格式(3.5)与隐格式(3.9a）是 $O\left(\tau +h^2\right)$ ，故从截断误差方面看，六点对称格式和Richardson格式有更大优越性.

(3) 稳定性. 在计算过程中, 由于初始数据有误差, 并且不可避免地有舍入误差, 因此人们自然关心这些误差传递下去, 是无限增长还是可以被控制? 这便是稳定性问题. 显然, 只有稳定的差分格式才是可用的.

作为例子, 我们考察Richardson 格式的稳定性. 由前面知道, Richardson 格式是显格式, 截断误差的阶是 $O\left(\tau^{2} + h^{2}\right)$ . 但从稳定性方面来看, 它是不可用的. 用 $e_{j}^{n}$ 表示 $u_{j}^{n}$ 的误差, 假定右端 $f_{j}^{n}$ 的计算是精确的, 则 $e_{j}^{n}$ 满足与 (3.16) 相应的齐方程:

$$
e _ {j} ^ {n + 1} = 2 r \left(e _ {j + 1} ^ {n} - 2 e _ {j} ^ {n} + e _ {j - 1} ^ {n}\right) + e _ {j} ^ {n - 1}. \tag {3.17}
$$

设误差只在初始层的原点 $(j = 0)$ 发生，即 $e_j^0 = \delta_{j0}\varepsilon (\varepsilon >0;\delta_{00} = 1,\delta_{j0} = 0$ 当 $j\neq 0$ 时)， $e_j^{-1} = 0$ ，而在以后计算中都是精确的，则初始误差的传播如表3.1.

表 3.1 $r = \frac{1}{2}$ 时 Richardson 格式的误差传播  

<table><tr><td rowspan="2">n</td><td colspan="9">j</td></tr><tr><td>-4</td><td>-3</td><td>-2</td><td>-1</td><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td></tr><tr><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>ε</td><td>0</td><td>0</td><td>0</td><td>0</td></tr><tr><td>1</td><td>0</td><td>0</td><td>0</td><td>ε</td><td>-2ε</td><td>ε</td><td>0</td><td>0</td><td>0</td></tr><tr><td>2</td><td>0</td><td>0</td><td>ε</td><td>-4ε</td><td>7ε</td><td>-4ε</td><td>ε</td><td>0</td><td>0</td></tr><tr><td>3</td><td>0</td><td>ε</td><td>-6ε</td><td>17ε</td><td>-24ε</td><td>17ε</td><td>-6ε</td><td>ε</td><td>0</td></tr><tr><td>4</td><td>ε</td><td>-8ε</td><td>31ε</td><td>-68ε</td><td>89ε</td><td>-68ε</td><td>31ε</td><td>-8ε</td><td>ε</td></tr><tr><td>5</td><td>-10ε</td><td>49ε</td><td>-144ε</td><td>273ε</td><td>-388ε</td><td>273ε</td><td>-144ε</td><td>49ε</td><td>-10ε</td></tr><tr><td>6</td><td>71ε</td><td>-260ε</td><td>641ε</td><td>-1096ε</td><td>1311ε</td><td>-1096ε</td><td>641ε</td><td>-260ε</td><td>71ε</td></tr></table>

从表中看出，误差随 $n\to \infty$ （ $h\rightarrow 0$ ）无限增长，所以差分格式不稳定.表中的计算虽然是就 $r = \frac{1}{2}$ 进行的，实际上对任何 $r > 0$ 都有类似现象，所以Richardson格式恒不稳定.

如果采用向前差分格式, 并取 $r = \frac{1}{2}$ , 则误差方程为

$$
e _ {j} ^ {n + 1} = \frac {1}{2} \left(e _ {j + 1} ^ {n} + e _ {j - 1} ^ {n}\right). \tag {3.18}
$$

此时误差逐渐衰减, 如表 3.2 所示. 显然如此的误差传递是允许的. 若限制 $0 < r \leqslant \frac{1}{2}$ , 则误差仍然衰减; 但当 $r > \frac{1}{2}$ 时, 误差也无限增长, 所以向前差分格式是条件稳定的.

表3.2 $r = \frac{1}{2}$ 时向前差分格式的误差传播  

<table><tr><td rowspan="2">n</td><td colspan="9">j</td></tr><tr><td>-4</td><td>-3</td><td>-2</td><td>-1</td><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td></tr><tr><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>ε</td><td>0</td><td>0</td><td>0</td><td>0</td></tr><tr><td>1</td><td>0</td><td>0</td><td>0</td><td>0.5ε</td><td>0</td><td>0.5ε</td><td>0</td><td>0</td><td>0</td></tr><tr><td>2</td><td>0</td><td>0</td><td>0.25ε</td><td>0</td><td>0.5ε</td><td>0</td><td>0.25ε</td><td>0</td><td>0</td></tr><tr><td>3</td><td>0</td><td>0.125ε</td><td>0</td><td>0.375ε</td><td>0</td><td>0.375ε</td><td>0</td><td>0.125ε</td><td>0</td></tr><tr><td>4</td><td>0.0625ε</td><td>0</td><td>0.25ε</td><td>0</td><td>0.375ε</td><td>0</td><td>0.25ε</td><td>0</td><td>0.0625ε</td></tr></table>

稳定性不仅对控制误差增长是重要的，而且也和收敛性有关，因此稳定性理论在数值求解非驻定问题中占有中心地位。

# 3.1.1 习题

1. (实习题) 就 $r = \frac{5}{11}, \frac{5}{9}$ 用显格式 3.7 作数值实验, 观察误差的增长规律, 并说明 $r = \frac{5}{11}$ 时稳定, $r = \frac{5}{9}$ 时不稳定.  
2. 将向前差分格式和向后差分格式作加权平均，得到下列格式：

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = \frac {a}{h ^ {2}} \left[ \theta \left(u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}\right) + (1 - \theta) \left(u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}\right) \right], \tag {3.19}
$$

其中 $0 \leqslant \theta \leqslant 1$ . 试计算其截断误差, 并证明当 $\theta = \frac{1}{2} - \frac{1}{12r}$ 时, 截断误差的阶最高 $(O(\tau^2 + h^4))$ .

3. 在Richardson格式(3.15)中以 $u_{j}^{n} = \frac{1}{2}\left(u_{j}^{n + 1} + u_{j}^{n - 1}\right)$ 代入左端，便得Dufort-Frankel格式：

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n - 1}}{2 \tau} = a \frac {u _ {j + 1} ^ {n} - u _ {j} ^ {n + 1} - u _ {j} ^ {n - 1} + u _ {j - 1} ^ {n}}{h ^ {2}}. \tag {3.20}
$$

试求其截断误差

4. 设有逼近热传导方程的带权三层差分格式：

$$
(1 + \theta) \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - \theta \frac {u _ {j} ^ {n} - u _ {j} ^ {n - 1}}{\tau} = a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}}, \tag {3.21}
$$

其中 $\theta \geqslant 0$ 试计算其截断误差，并证明当 $\theta = \frac{1}{2} +\frac{1}{12r}$ 时，截断误差的阶最高 $(O(\tau^2 +h^4))$

# 3.2 稳定性与收敛性

# 3.2.1 稳定性概念

前节引进的二层差分格式, 均可用矩阵和向量的记号表成

$$
\boldsymbol {A} \boldsymbol {U} ^ {n + 1} = \boldsymbol {B} \boldsymbol {U} ^ {n} + \tau \boldsymbol {F}, \tag {3.22}
$$

其中 $\pmb{U}^{n} = \left(u_{1}^{n},u_{2}^{n},\dots ,u_{J - 1}^{n}\right)^{\mathrm{T}},\pmb {F} = \left(f_{1},f_{2},\dots ,f_{J - 1}\right)^{\mathrm{T}},\pmb{A}$ 和 $\pmb{B}$ 是 $(J - 1)\times (J - 1)$ 矩阵.假定 $\mathbf{A}$ 有逆，并令

$$
C = A ^ {- 1} B, \tag {3.23}
$$

则可将(3.22)化为

$$
\boldsymbol {U} ^ {n + 1} = \boldsymbol {C} \boldsymbol {U} ^ {n} + \tau \boldsymbol {A} ^ {- 1} \boldsymbol {F}, \tag {3.24}
$$

其中 $C$ 称为增长矩阵

例如对于向前差分格式, $A = I((J - 1)$ 阶单位矩阵), $B = (1 - 2r)I + rS,$ 其中

$$
S = \left[ \begin{array}{c c c c c c} 0 & 1 & & & & \\ 1 & 0 & & & & \\ & & \ddots & & \\ & & \ddots & & \\ & & & & 1 \\ & & & & 1 & 0 \end{array} \right] _ {(J - 1) \times (J - 1)} \tag {3.25}
$$

故 $C = (1 - 2r)I + rS.$ 对于向后差分格式， $\pmb {A} = (1 + 2r)\pmb {I} - r\pmb {S},\pmb {B} = \pmb{I}$ ，故 $C =$ $[(1 + 2r)\pmb {I} - r\pmb {S}]^{-1}$ .对于六点对称格式， $\pmb {A} = (1 + r)\pmb {I} - \frac{r}{2}\pmb {S},\pmb {B} = (1 - r)\pmb {I} + \frac{r}{2}\pmb{S},$ 故 $C = \left[(1 + r)\pmb {I} - \frac{r}{2}\pmb {S}\right]^{-1}\left[(1 - r) + \frac{r}{2}\pmb {S}\right].$

至于一般的三层或多层格式，总可适当引进新变量化成二层格式。例如Richardson格式，其矩阵形式为

$$
\boldsymbol {U} ^ {n + 1} = 2 r (\boldsymbol {S} - 2 \boldsymbol {I}) \boldsymbol {U} ^ {n} + \boldsymbol {U} ^ {n - 1}. \tag {3.26}
$$

令 $\pmb{W}^{n} = \left(\pmb{U}^{n},\pmb{U}^{n - 1}\right)^{\mathrm{T}}$ ，则化为

$$
\boldsymbol {W} ^ {n + 1} = \boldsymbol {C W} ^ {n}, \tag {3.27}
$$

其中

$$
\boldsymbol {C} = \left[ \begin{array}{c c} 2 r (\boldsymbol {S} - 2 \boldsymbol {I}) & \boldsymbol {I} \\ \boldsymbol {I} & \boldsymbol {O} \end{array} \right]. \tag {3.28}
$$

我们仅讨论系数及右端与时间 $t$ 无关的线性抛物型方程，所以(3.22)中的 $\mathbf{A},\mathbf{B}$ 和 $\pmb{F}$ 均不依赖 $n,$ 但 $A,B$ 可依赖步长 $h$ 和 $\tau .$ 我们要求 $h,\tau$ 之间满足一定关系，设为 $h = g(\tau)$ ，其中 $g(\tau)$ 连续且满足 $g(0) = 0.$ 于是 $\mathbf{A} = \mathbf{A}(\tau),\mathbf{B} = \mathbf{B}(\tau),\mathbf{C} = \mathbf{C}(\tau).$

先讨论按初值稳定.此时 $F = 0$

$$
\boldsymbol {U} ^ {n + 1} = \boldsymbol {C} (\tau) \boldsymbol {U} ^ {n} = \dots = [ \boldsymbol {C} (\tau) ] ^ {n + 1} \boldsymbol {U} ^ {0}. \tag {3.29}
$$

我们说差分格式 (3.22) 按初值稳定, 如果存在 $\tau_0 > 0$ 和常数 $K > 0$ , 使不等式

$$
\left\| \boldsymbol {U} ^ {n + 1} \right\| = \left\| [ \boldsymbol {C} (\tau) ] ^ {n + 1} \boldsymbol {U} ^ {0} \right\| \leqslant K \left\| \boldsymbol {U} ^ {0} \right\| \tag {3.30}
$$

对一切 $U^0\in \mathbb{R}^{J - 1},0 <   \tau \leqslant \tau_0$ 和 $0 <   n\tau \leqslant T$ 成立.这里 $\| \cdot \|$ 是 $\mathbb{R}^{J - 1}$ 中的某一种范数，一般取

$$
\| \boldsymbol {U} \| ^ {2} = \sum_ {j = 1} ^ {J - 1} h u _ {j} ^ {2}.
$$

显然差分格式 (3.22) 按初值稳定, 当且只当

$$
\| \boldsymbol {C} ^ {n} (\tau) \| \leqslant K, \quad 0 <   \tau \leqslant \tau_ {0}, \quad 0 <   n \tau \leqslant T. \tag {3.31}
$$

其次讨论按右端稳定. 此时认为初值没有误差, 即 $U^0 = 0$ . 我们说差分格式 (3.22) 按右端稳定, 如果存在 $\tau_0 > 0$ 和常数 $K > 0$ , 使不等式

$$
\left\| U ^ {n + 1} \right\| \leqslant K \| F \|
$$

对一切 $0 < \tau \leqslant \tau_0$ 和 $0 < n\tau \leqslant T$ 成立, 其中 $\pmb{U}^{n}$ 是下列方程的解:

$$
\boldsymbol {U} ^ {n + 1} = \boldsymbol {C} (\tau) \boldsymbol {U} ^ {n} + \tau \boldsymbol {A} ^ {- 1} \boldsymbol {F}, \quad \boldsymbol {U} ^ {0} = \mathbf {0}. \tag {3.32}
$$

反复利用递推式 (3.32), 得

$$
\begin{array}{l} \boldsymbol {U} ^ {n + 1} = \boldsymbol {C} (\tau) \boldsymbol {U} ^ {n} + \tau \boldsymbol {A} ^ {- 1} \boldsymbol {F} \\ = C (\tau) \left[ C (\tau) U ^ {n - 1} + \tau A ^ {- 1} F \right] + \tau A ^ {- 1} F \\ = C ^ {2} (\tau) U ^ {n - 1} + \tau C (\tau) A ^ {- 1} F + \tau A ^ {- 1} F \\ = C ^ {2} (\tau) \left[ C (\tau) U ^ {n - 2} + \tau A ^ {- 1} F \right] + \tau \left[ C (\tau) + I \right] A ^ {- 1} F \\ \begin{array}{c c c} \cdot & \cdot & \cdot \\ \cdot & \cdot & \cdot \end{array} \\ = \tau \left[ \boldsymbol {C} ^ {n} (\tau) + \boldsymbol {C} ^ {n - 1} (\tau) + \dots + \boldsymbol {C} (\tau) + \boldsymbol {I} \right] \boldsymbol {A} ^ {- 1} \boldsymbol {F}. \\ \end{array}
$$

设 $\| A^{-1}\| \leqslant K^{\prime}$ ，又差分格式按初值稳定，即存在常数 $K''$ 使 $\| C^n (\tau)\| \leqslant K''$ ，则

$$
\left\| U ^ {n + 1} \right\| \leqslant \tau (n + 1) K ^ {\prime} K ^ {\prime \prime} \| F \| \leqslant T K ^ {\prime} K ^ {\prime \prime} \| F \|.
$$

取 $K = TK^{\prime}K^{\prime \prime}$ ，即知格式按右端稳定

如果右端与时间有关，即 $F = F^n$ ，以上推导仍成立，只需用 $\sup_n\| F^n\|$ 代替上述不等式右端的 $\| F\|$ 总之，若 $\left\| A^{-1}(\tau)\right\| \leqslant K^{\prime}$ ，则由格式按初值稳定可推出它按右端稳定为检验格式按初值稳定，需检验不等式(3.31)，即矩阵族

$$
\left\{\boldsymbol {C} ^ {n} (\tau): 0 <   \tau \leqslant \tau_ {0}, 0 <   n \tau \leqslant T \right\} \tag {3.33}
$$

一致有界. 往后, 我们所述的稳定均指按初值稳定

# 3.2.2 判别稳定性的直接估计法（矩阵法）

判别矩阵族 (3.33) 的一致有界性是一个困难问题, 只在某些特殊情形才能给出解答.

命题3.1（必要条件）以 $\rho (C)$ 表示矩阵 $C(\tau)$ 的谱半径，则差分格式稳定的必要条件是存在与 $\tau$ 无关的常数 $M$ 使

$$
\rho (\boldsymbol {C}) \leqslant 1 + M \tau \quad (\rho (\boldsymbol {C}) \leqslant 1 + O (\tau)). \tag {3.34}
$$

证明 由(3.31)知，

$$
\rho^ {n} (\boldsymbol {C}) \leqslant \| \boldsymbol {C} ^ {n} \| \leqslant K, \quad 0 <   n \leqslant \frac {T}{\tau}, \quad 0 <   \tau \leqslant \tau_ {0}.
$$

不妨设 $K > 1$ ，并取 $n = \left[\frac{T}{\tau}\right]\left(\left[\frac{T}{\tau}\right]$ 表示 $\frac{T}{\tau}$ 的整数部分），则

$$
\rho (C) \leqslant K ^ {\frac {1}{n}} \leqslant K ^ {\frac {\tau}{T - \tau}} = \mathrm {e} ^ {\frac {\tau}{T - \tau} \ln K} \leqslant \mathrm {e} ^ {\frac {\ln K}{T - \tau_ {0}} \tau} = 1 + O (\tau).
$$

命题3.2（充分条件）若 $C(\tau)$ 是正规矩阵，即 $C$ 和它的共轭转置 $C^*$ 乘积可交换： $CC^{*} = C^{*}C$ ，则(3.34）也是差分格式稳定的充分条件。

证明 因为此时 $\|C(\tau)\| = \rho(C)$ , 由 (3.34),

$$
\begin{array}{l} \| \boldsymbol {C} ^ {n} (\tau) \| \leqslant \| \boldsymbol {C} (\tau) \| ^ {n} = \rho^ {n} (\boldsymbol {C}) \leqslant (1 + M \tau) ^ {n} \\ \leqslant (1 + M \tau) ^ {\frac {T}{\tau}} \leqslant K <   \infty . \\ \end{array}
$$

推论3.1 若 $S$ 是对称矩阵， $C(\tau)$ 是矩阵 $S$ 的实系数有理函数： $C(\tau) = R(S)$ ，则差分格式稳定的充要条件是

$$
\max  _ {j} \left| R \left(\lambda_ {j} ^ {S}\right) \right| \leqslant 1 + M \tau ,
$$

其中 $\lambda_j^S$ 是 $s$ 的特征值

注意 $R(S)$ 是实数和矩阵 $\pmb{S}$ 的四则运算，因此矩阵 $C(\tau)$ 也对称，其特征值是 $R(\lambda_j^S)$ 特别地，当 $\pmb{S}$ 是形如（3.25）的矩阵时，其特征值

$$
\lambda_ {j} ^ {S} = 2 \cos j \pi h, \quad j = 1, 2, \dots , J - 1, \quad h = \frac {l}{J},
$$

特征向量 $U^{j}$ 的分量 $u_{k}^{j} = \sin jk\pi h,k = 1,2,\dots ,J - 1.$

例3.1 对向前差分格式（以下设(3.4)中的 $l = 1$ ）， $C = (1 - 2r)\pmb{I} + r\pmb{S}$ ，

$$
\lambda_ {j} ^ {C} = 1 - 2 r + 2 r \cos j \pi h = 1 - 4 r \sin^ {2} \frac {j \pi h}{2}.
$$

为使 $|\lambda_j^C |\leqslant 1 + M\tau$ 或

$$
- 1 - M \tau \leqslant \lambda_ {j} ^ {C} = 1 - 4 r \sin^ {2} \frac {j \pi h}{2} \leqslant 1 + M \tau ,
$$

必须且只需

$$
4 r \sin^ {2} \frac {j \pi h}{2} \leqslant 2 + M \tau , \quad j = 1, 2, \dots , J - 1,
$$

从而 $4r \leqslant 2, r \leqslant \frac{1}{2}$ . 所以向前差分格式当 $r \leqslant \frac{1}{2}$ 时稳定, 当 $r > \frac{1}{2}$ 时不稳定.

例3.2 对向后差分格式, $C = [(1 + 2r)\pmb {I} - r\pmb{S}]^{-1}$

$$
\lambda_ {j} ^ {\mathbf {C}} = \left[ (1 + 2 r) - 2 r \cos j \pi h \right] ^ {- 1} = \left[ 1 + 2 r (1 - \cos j \pi h) \right] ^ {- 1} \leqslant 1,
$$

故对任何 $r > 0$ 稳定，即恒稳定或绝对稳定

例3.3 对六点对称格式，

$$
\boldsymbol {C} = \left[ (1 + r) \boldsymbol {I} - \frac {r}{2} \boldsymbol {S} \right] ^ {- 1} \left[ (1 - r) \boldsymbol {I} + \frac {r}{2} \boldsymbol {S} \right],
$$

$$
\lambda_ {j} ^ {C} = \frac {1 - 2 r \sin^ {2} \frac {j \pi h}{2}}{1 + 2 r \sin^ {2} \frac {j \pi h}{2}}, \quad j = 1, 2, \dots , J - 1,
$$

故对任何 $r > 0$ 有 $\left|\lambda_j^{\mathcal{C}}\right| \leqslant 1$ ，因此六点对称格式恒稳定

例3.4 对Richardson格式，

$$
C = \left[ \begin{array}{c c} 2 r (S - 2 I) & I \\ I & O \end{array} \right]
$$

是对称矩阵.设 $\lambda$ 为 $C$ 的特征值， $\pmb {W} = (\omega_{1},\omega_{2})^{\mathrm{T}}$ 为相应的特征向量，即 $CW = \lambda W,$ 或

$$
2 r (S - 2 I) \omega_ {1} + \omega_ {2} = \lambda \omega_ {1}, \quad \omega_ {1} = \lambda \omega_ {2}.
$$

显然 $\omega_{2} \neq 0$ . 利用第二方程消去 $\omega_{1}$ , 得

$$
2 \lambda r (S - 2 I) \omega_ {2} + \omega_ {2} = \lambda^ {2} \omega_ {2},
$$

从而

$$
S \omega_ {2} = \left(2 + \frac {\lambda}{2 r} - \frac {1}{2 \lambda r}\right) \omega_ {2}.
$$

可见 $\mu = 2 + \frac{\lambda}{2r} -\frac{1}{2\lambda r}$ 是 $s$ 的特征值.于是 $\lambda$ 满足方程

$$
\lambda^ {2} + 2 r (2 - \mu) \lambda - 1 = 0 (\mu = 2 \cos j \pi h),
$$

或

$$
\lambda^ {2} + \left(8 r \sin^ {2} \frac {j \pi h}{2}\right) \lambda - 1 = 0.
$$

其根的按模最大值

$$
\max  _ {j} \left\{\left| \lambda_ {1} ^ {j} \right|, \left| \lambda_ {2} ^ {j} \right| \right\} = \max  _ {j} \left\{\left| 4 r \sin^ {2} \frac {j \pi h}{2} + \sqrt {1 6 r ^ {2} \sin^ {4} \frac {j \pi h}{2} + 1} \right| \right\}
$$

$>r+\sqrt{1+r^{2}}>1+r,$ 对任意 $r>0$

所以Richardson格式恒不稳定

例3.5 带第二边值条件的向前差分格式：

$$
u _ {j} ^ {n + 1} = r u _ {j + 1} ^ {n} + (1 - 2 r) u _ {j} ^ {n} + u _ {j - 1} ^ {n}, \quad j = 1, 2, \dots , J - 1, \tag {3.35a}
$$

$$
u _ {0} ^ {n} = 0, \quad u _ {J} ^ {n} = u _ {J - 1} ^ {n}, \tag {3.35b}
$$

其中 $u_{J}^{n} = u_{J - 1}^{n}$ 是第二边值条件 $u_{x}(1,t) = 0$ 的差分近似.利用(3.35b）消去（3.35a）中的 $u_0^n,u_J^n$ ，则知增长矩阵

$$
C = \left[ \begin{array}{c c c c c} 1 - 2 r & r & & & \\ r & 1 - 2 r & r & & \\ & \ddots & \ddots & \ddots & \\ & & r & 1 - 2 r & r \\ C (r) & & & r & 1 - r \end{array} \right] = I - r D,
$$

其中

$$
\boldsymbol {D} = \left[ \begin{array}{c c c c c} 2 & - 1 & & & \\ - 1 & 2 & - 1 & & \\ & \ddots & \ddots & \ddots & \\ & & - 1 & 2 & - 1 \\ & & & - 1 & 1 \end{array} \right]
$$

是对称矩阵，其特征值 $\lambda$ 是实数.求 $D$ 的特征值等同于解下列差分算子的特征问题：

$$
l _ {h} u _ {j} = - \left(u _ {j + 1} + u _ {j - 1}\right) = (\lambda - 2) u _ {j}, \quad j = 1, 2, \dots , J - 1, \quad u _ {0} = 0, \quad u _ {J} = u _ {J - 1}.
$$

以 $u_{j} = z^{j}$ 代入上式，得 $z^2 +(\lambda -2)z + 1 = 0$ ，其解

$$
z _ {1}, z _ {2} = \frac {1}{2} \left(2 - \lambda \pm \sqrt {\lambda^ {2} - 4 \lambda}\right).
$$

为使解满足左边值条件，它的二根 $z_{1},z_{2}$ 不能是实的，而是一对共轭复根 $z_{1} = \bar{z}_{2} = \mathrm{e}^{\mathrm{i}\theta}$ 且 $\cos \theta = 1 - \frac{\lambda}{2},\sin \theta = \frac{\sqrt{4\lambda - \lambda^2}}{2}$ 一般解为 $u_{j} = c_{1}\cos j\theta +c_{2}\sin j\theta .$ 由边值条件 $u_0 = 0,u_J = u_{J - 1}$ 知 $c_{1} = 0,\sin J\theta = \sin (J - 1)\theta ,$ 从而 $\theta = \frac{(2j - 1)\pi}{2J - 1},D$ 的特征值为 $\lambda = 2 - 2\cos \theta$ ，即

$$
\begin{array}{l} \lambda_ {j} = 2 - 2 \cos \frac {(2 j - 1) \pi}{2 J - 1}, \\ = 4 \sin^ {2} \frac {(2 j - 1) \pi}{2 (2 J - 1) \pi}, \quad j = 1, 2, \dots , J - 1. \\ \end{array}
$$

所以 $C$ 的特征值

$$
\lambda_ {j} ^ {C} = 1 - 4 r \sin^ {2} \frac {2 j - 1}{2 (2 J - 1)} \pi , \quad j = 1, 2, \dots , J - 1.
$$

为使 $|\lambda_j^C |\leqslant 1$ 必须且只需 $r\leqslant \frac{1}{2}$ 又因 $c$ 对称，故格式（3.35a）（3.35b）的稳定条件是 $r\leqslant \frac{1}{2}.$ 四i

# 3.2.3 收敛性与敛速估计

考虑热传导方程的初边值问题：

$$
\left\{ \begin{array}{l l} L u = \frac {\partial u}{\partial t} - a \frac {\partial^ {2} u}{\partial x ^ {2}} = f (x), & 0 <   x <   l, \quad 0 <   t \leqslant T, \\ u (x, 0) = \varphi (x), & u (0, t) = u (l, t) = 0. \end{array} \right. \tag {3.36}
$$

相应的差分格式为

$$
\left\{ \begin{array}{l l} L _ {h} u _ {j} ^ {n} = f _ {j}, & j = 1, 2, \dots , J - 1, \quad n = 0, 1, \dots , N - 1, \\ u _ {j} ^ {0} = \varphi_ {j}, & u _ {0} ^ {n} = u _ {J} ^ {n} = 0. \end{array} \right. \tag {3.37}
$$

其向量形式如 (3.24).

差分逼近的截断误差

$$
R _ {j} ^ {n} (u) = L _ {h} u \left(x _ {j}, t _ {n}\right) - [ L u ] _ {j} ^ {n}, \tag {3.38}
$$

$u(x,t)$ 是 $0\leqslant x\leqslant l,0\leqslant t\leqslant T$ 上的任一充分光滑函数.称差分算子 $L_{h}$ 是边值问题(3.36）的相容逼近，如果相容条件

$$
\lim  _ {\tau \rightarrow 0} \| \boldsymbol {R} ^ {n} \| = 0 \quad (\| \boldsymbol {R} ^ {n} \| = o (1)) \tag {3.39}
$$

成立, 其中 $R^n$ 是分量为 $R_j^n(u)$ 的向量, $\|\cdot\|$ 是 $\mathbb{R}^{J-1}$ 中的范数 (参看 (2.11)-(2.14)).

先对差分解作出某种估计. 将 (3.37) 的解分解为 $u_{j}^{n} = v_{j}^{n} + w_{j}^{n}$ , 其中 $v_{j}^{n}$ 满足零初值和非齐右端方程:

$$
\boldsymbol {V} ^ {n + 1} = \boldsymbol {C} (\tau) \boldsymbol {V} ^ {n} + \tau \boldsymbol {A} ^ {- 1} \boldsymbol {F} \quad \left(\boldsymbol {V} ^ {0} = \boldsymbol {0}\right),
$$

而 $w_{j}^{n}$ 满足非零初值和齐右端方程：

$$
\boldsymbol {W} ^ {n + 1} = \boldsymbol {C} (\tau) \boldsymbol {W} ^ {n} \quad \left(\boldsymbol {W} ^ {0} = \boldsymbol {U} ^ {0}\right),
$$

其中 $V^n, W^n$ 依次为以 $v_j^n, w_j^n$ 为分量的向量. 若差分格式按初值稳定, 则亦按右端稳定, 于是有常数 $K_1$ 和 $K_2$ , 使

$$
\| \boldsymbol {V} ^ {n} \| \leqslant K _ {1} \| \boldsymbol {F} \|, \| \boldsymbol {W} ^ {n} \| \leqslant K _ {2} \| \boldsymbol {W} ^ {0} \| = K _ {2} \| \boldsymbol {U} ^ {0} \|.
$$

这样

$$
\left\| \boldsymbol {U} ^ {n} \right\| \leqslant K \left(\left\| \boldsymbol {U} ^ {0} \right\| + \left\| \boldsymbol {F} \right\|\right), \quad K = \max  \left\{K _ {1}, K _ {2} \right\}. \tag {3.40}
$$

现在估计差分解的误差. 设 $u(x, t)$ 是热传导方程 (3.36) 的解, $u_{j}^{n}$ 是差分方程 (3.37) 的解. 误差 $e_{j}^{n} = u(x_{j}, t_{n}) - u_{j}^{n} = [u]_{j}^{n} - u_{j}^{n}$ . 我们有

$$
\begin{array}{l} R _ {j} ^ {n} (u) = L _ {h} [ u ] _ {j} ^ {n} - \left[ L u \right] _ {j} ^ {n} = L _ {h} [ u ] _ {j} ^ {n} - f _ {j} \\ = L _ {h} [ u ] _ {j} ^ {n} - L _ {h} u _ {j} ^ {n} = L _ {h} e _ {j} ^ {n}, \\ \end{array}
$$

即误差 $e_j^n$ 满足差分方程：

$$
L _ {h} e _ {j} ^ {n} = R _ {j} ^ {n} (u), \quad e _ {j} ^ {0} = 0, \quad j = 1, 2, \dots , J - 1.
$$

由(3.24)知其向量形式为

$$
\boldsymbol {E} ^ {n + 1} = \boldsymbol {C} (\tau) \boldsymbol {E} ^ {n} + \tau \boldsymbol {A} ^ {- 1} \boldsymbol {R} ^ {n},
$$

这里 $E^n, R^n$ 依次为以 $e_j^n, R_j^n$ 为分量的向量. 由估计式 (3.40) 得

$$
\| \boldsymbol {E} ^ {n} \| \leqslant K \sup  _ {n} \| \boldsymbol {R} ^ {n} \|. \tag {3.41}
$$

若相容条件 (3.39) 成立, 则

$$
\lim  _ {\tau \rightarrow 0} \| \boldsymbol {E} ^ {n} \| = \lim  _ {\tau \rightarrow 0} \| \boldsymbol {u} ^ {n} - \boldsymbol {U} ^ {n} \| = 0,
$$

其中 $\pmb{u}^n$ 表示以 $u(x_{j},t_{n})$ 为分量的向量.这证明了如下

定理3.1 若差分方程满足相容条件且按初值稳定，则差分解收敛到热传导方程的解且有误差估计式(3.41).

推论3.2 当网比 $r \leqslant \frac{1}{2}$ 时，向前差分格式的解有收敛阶 $O\left(\tau + h^{2}\right)$ . 对任何网比 $r > 0$ ，向后差分格式的解有收敛阶 $O\left(\tau + h^{2}\right)$ ，六点对称格式的解有收敛阶 $O\left(\tau^{2} + h^{2}\right)$ .

注3.1实际上，定理3.1及其证明对更一般的非驻定偏微分方程和差分格式也成立（参见[28]的第3章）

# 3.2.4 习题

1. 求证差分格式 (3.19) 当 $\frac{1}{2} \leqslant \theta \leqslant 1$ 时恒稳定, 当 $0 \leqslant \theta < \frac{1}{2}$ 时稳定的充要条件是

$$
r \leqslant \frac {1 - 2 \theta}{2}.
$$

2. 证明如下格式恒稳定：

$$
\begin{array}{l} \frac {1}{1 2} \frac {u _ {j + 1} ^ {n + 1} - u _ {j + 1} ^ {n}}{\tau} + \frac {5}{6} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + \frac {1}{1 2} \frac {u _ {j - 1} ^ {n + 1} - u _ {j - 1} ^ {n}}{\tau} \\ = a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1} + u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{2 h ^ {2}} (a > 0). \tag {3.42} \\ \end{array}
$$

# 3.3 Fourier 方法

前节介绍的判别稳定性的直接估计法，原则上可用于一般非驻定问题，但只在某些简单情形才能估计矩阵族 $C^n(\tau)$ 的范数，这里遇到的主要困难之一是矩阵 $C(\tau)$ 的阶 $(J - 1)$ 随 $\tau \to 0$ 而无限增大。本节仅限于讨论常系数线性非驻定方程的纯初值问题和带周期边值条件的混合问题，此时可用Fourier方法（Fourier积分和Fourier级数）将空间变量和时间变量分离，从而将差分方程的稳定性归结为有限阶矩阵族的一致有界性。

先考虑线性常系数一维抛物型方程，具初值和周期（设周期为 $l$ ）边值条件。逼近它的二层差分方程的一般形式为

$$
\sum_ {m \in \mathcal {N} _ {1}} a _ {m} u _ {j + m} ^ {n + 1} = \sum_ {m \in \mathcal {N} _ {0}} b _ {m} u _ {j + m} ^ {n}, \quad j = 0, 1, \dots , J - 1 \tag {3.43}
$$

(只考虑按初值稳定，故可设非齐项等于零).这是在空间网点 $x_{j}$ 处的差分方程， $\mathcal{N}_0$ 和 $\mathcal{N}_{1}$ 是包含0及其附近的正负整数的有限集合， $a_{m}$ 和 $b_{m}$ 不依赖 $j$ 但可能和 $\tau$ 有关.例如对向前差分格式

$$
\begin{array}{l} u _ {j} ^ {n + 1} = r u _ {j + 1} ^ {n} + (1 - 2 r) u _ {j} ^ {n} + r u _ {j - 1} ^ {n}, \\ \mathcal {N} _ {0} = \{- 1, 0, 1 \}, \quad \mathcal {N} _ {1} = \{0 \}, \quad b _ {- 1} = b _ {1} = r, \quad b _ {0} = 1 - 2 r, \quad a _ {0} = 1. \\ \end{array}
$$

对向后差分格式

$$
\begin{array}{l} - r u _ {j + 1} ^ {n + 1} + (1 + 2 r) u _ {j} ^ {n + 1} - r u _ {j - 1} ^ {n + 1} = u _ {j} ^ {n}, \\ \mathcal {N} _ {0} = \{0 \}, \quad \mathcal {N} _ {1} = \{- 1, 0, 1 \}, \quad b _ {0} = 1, \quad a _ {- 1} = a _ {1} = - r, \quad a _ {0} = 1 + 2 r. \\ \end{array}
$$

对六点对称格式

$$
\begin{array}{l} - \frac {r}{2} u _ {j + 1} ^ {n + 1} + (1 + r) u _ {j} ^ {n + 1} - \frac {r}{2} u _ {j - 1} ^ {n + 1} = \frac {r}{2} u _ {j + 1} ^ {n} + (1 - r) u _ {j} ^ {n} + \frac {r}{2} u _ {j - 1} ^ {n}, \\ \mathcal {N} _ {0} = \mathcal {N} _ {1} = \{- 1, 0, 1 \}, \quad b _ {- 1} = b _ {1} = \frac {r}{2}, \quad b _ {0} = 1 - r, \quad a _ {- 1} = a _ {1} = - \frac {r}{2}, \quad a _ {0} = 1 + r. \\ \end{array}
$$

由于是周期边值条件 $(u_0^n = u_j^n)$ , 故可将 $u_{j}^{n}$ 周期开拓使其对一切 $j = 0, \pm 1, \dots$ 有意义, 且方程 (3.43) 对所有整数 $j$ 成立. 为了应用 Fourier 方法, 我们再将 $u_{j}^{n} = u^{n}(x_{j})$ 开拓为 $(- \infty, \infty)$ 上的 $u^{n}(x)$ . 为此, 取半整数点 $x_{j + \frac{1}{2}} = x_{j} + \frac{1}{2} h, j = 0, \pm 1, \dots$ , 并用如下阶梯函数逼近初始函数 $\varphi(x)$ :

$$
u ^ {0} (x) = \varphi (x _ {j}), \quad x _ {j - \frac {1}{2}} <   x <   x _ {j + \frac {1}{2}},
$$

其中 $j = 0,\pm 1,\dots$ 再将(3.43)看成在任一 $x_{j} = x\in (-\infty ,\infty)$ 成立，则得具连续变量的差分解 $u^{n}(x)$ .显然 $u^n (x)$ 仍是 $\mathcal{X}$ 的周期函数（周期为 $l$ ),且

$$
u ^ {n} (x) = u _ {j} ^ {n}, \quad x _ {j - \frac {1}{2}} <   x <   x _ {j + \frac {1}{2}},
$$

其中 $j = 0,\pm 1,\dots$ 显然 $u^{n}(x)$ 于 $(0,l)$ 平方可积，因此属于空间 $L^2 (0,l)$ ，其范数

$$
\| u ^ {n} \| _ {L ^ {2}} ^ {2} = \int_ {0} ^ {l} | u ^ {n} (x) | ^ {2} \mathrm {d} x.
$$

此外还有

$$
\begin{array}{l} \| U ^ {n} \| ^ {2} = \sum_ {j = 0} ^ {J - 1} h \left(u _ {j} ^ {n}\right) ^ {2} = \sum_ {j = 1} ^ {J - 1} h \left(u _ {j} ^ {n}\right) ^ {2} + \frac {h}{2} \left\{\left(u _ {0} ^ {n}\right) ^ {2} + \left(u _ {N} ^ {n}\right) ^ {2} \right\} \\ = \int_ {0} ^ {l} | u ^ {n} (x) | ^ {2} d x = \| u ^ {n} (x) \| _ {L ^ {2}} ^ {2}. \\ \end{array}
$$

这样，我们就可将Fourier方法用于具连续空间变量的差分方程

$$
\sum_ {m \in \mathcal {N} _ {1}} a _ {m} u ^ {n + 1} (x + x _ {m}) = \sum_ {m \in \mathcal {N} _ {0}} b _ {m} u ^ {n} (x + x _ {m}). \tag {3.44}
$$

将 $u^n (x)$ 展成Fourier级数：

$$
u ^ {n} (x) = \sum_ {p = - \infty} ^ {\infty} v _ {p} ^ {n} \exp \left(\mathrm {i} \frac {2 p \pi}{l} x\right), \tag {3.45}
$$

#

$$
v _ {p} ^ {n} = \frac {1}{l} \int_ {0} ^ {l} u ^ {n} (x) \exp \left(- \mathrm {i} \frac {2 p \pi}{l} x\right) \mathrm {d} x, \quad p = 0, \pm 1, \dots . \tag {3.46}
$$

人物简介我们有 Parseval 等式：

$$
\left\| u ^ {n} (x) \right\| _ {L ^ {2}} ^ {2} = l \sum_ {p = - \infty} ^ {\infty} \left| v _ {p} ^ {n} \right| ^ {2}. \tag {3.47}
$$

把(3.45)代到(3.44)，得

$$
\begin{array}{l} \sum_ {p = - \infty} ^ {\infty} v _ {p} ^ {n + 1} \left[ \sum_ {m \in \mathcal {N} _ {1}} a _ {m} \exp \left(\mathrm {i} \frac {2 p \pi}{l} x _ {m}\right) \right] \exp \left(\mathrm {i} \frac {2 p \pi}{l} x\right) \tag {3.48} \\ = \sum_ {p = - \infty} ^ {\infty} v _ {p} ^ {n} \left[ \sum_ {m \in \mathcal {N} _ {0}} b _ {m} \exp \left(\mathrm {i} \frac {2 p \pi}{l} x _ {m}\right) \right] \exp \left(\mathrm {i} \frac {2 p \pi}{l} x\right). \\ \end{array}
$$

比较对应项的系数，得

$$
v _ {p} ^ {n + 1} = G (p h, \tau) v _ {p} ^ {n}, \tag {3.49}
$$

其中

$$
G (p h, \tau) = \left[ \sum_ {m \in \mathcal {N} _ {1}} a _ {m} \exp \left(\mathrm {i} \frac {2 m \pi}{l} p h\right) \right] ^ {- 1} \left[ \sum_ {m \in \mathcal {N} _ {0}} b _ {m} \exp \left(\mathrm {i} \frac {2 m \pi}{l} p h\right) \right]. \tag {3.50}
$$

将(3.49)代到(3.47)，则

$$
\left\| u ^ {n} (x) \right\| _ {L ^ {2}} ^ {2} = l \sum_ {p = - \infty} ^ {\infty} \left| G (p h, \tau) v _ {p} ^ {n - 1} \right| ^ {2} = l \sum_ {p = - \infty} ^ {\infty} \left| G ^ {n} (p h, \tau) v _ {p} ^ {0} \right| ^ {2}. \tag {3.51}
$$

若差分格式稳定，则有常数 $K > 0$ 使

$$
\left\| u ^ {n} (x) \right\| _ {L ^ {2}} ^ {2} = l \sum_ {p = - \infty} ^ {\infty} \left| G ^ {n} (p h, \tau) v _ {p} ^ {0} \right| ^ {2} \leqslant l K ^ {2} \left\| u ^ {0} (x) \right\| _ {L ^ {2}} ^ {2}.
$$

由于阶梯函数类 $\{u^0 (x)\}$ 于 $L^2 (0,l)$ 稠密，取 $u^0 (x)$ 使 $\| u^0 (x)\|_{L^2} = 1$ ，其Fourier系数 $v_{p}^{0} = 1,v_{q}^{0} = 0(q\neq p),$ 则由上式得

$$
\left| G ^ {n} \left(p h, \tau\right) \right| \leqslant K, \quad 0 <   \tau \leqslant \tau_ {0}, \quad 0 <   n \tau \leqslant T, \tag {3.52}
$$

即 $G^{n}(ph,\tau)$ 一致有界.反之,若 $G^{n}(ph,\tau)$ 一致有界，则由(3.51）得

$$
\left\| u ^ {n} (x) \right\| _ {L ^ {2}} ^ {2} \leqslant K ^ {2} l \sum_ {p = - \infty} ^ {\infty} \left| v _ {p} ^ {0} \right| ^ {2} = K ^ {2} \left\| u ^ {0} (x) \right\| _ {L ^ {2}} ^ {2}, \tag {3.53}
$$

从而差分格式按初值稳定

往后我们称 $G(ph,\tau)$ 为增长因子（amplification factor).不妨设 $K\geqslant 1$ 显然不等式(3.52）又等价于 $\left|G(ph,\tau)\right|\leqslant K^{\frac{1}{n}}$ ，取 $n = \left[\frac{T}{\tau}\right]$ (取整数部分)，则 $\frac{1}{n}\leqslant \frac{\tau}{T - \tau}$ ，所以 $|G(ph,\tau)|\leqslant K^{\frac{\tau}{T - \tau}}\leqslant \exp \left(\frac{\tau}{T - \tau}\ln K\right)$ ，于是

$$
| G (p h, \tau) | \leqslant 1 + M \tau . \tag {3.54}
$$

上式也称为 von Neumann 条件. 综合上述, 我们得

命题3.3 差分格式(3.43)稳定 $\Longleftrightarrow G^{n}(ph,\tau)$ 一致有界 $\Longleftrightarrow$ von Neumann条件(3.54)成立.

注3.2 注意(3.50)中 $ph = x_{p}$ 是空间网点, $G(x_{p},\tau)$ 关于 $x_{p},\tau$ 连续, 关于 $x_{p}$ 是以 $l$ 为周期的函数, 所以只需就 $p = 0,1,\dots ,J - 1$ 研究 $G^{n}(x_{p},\tau)$ 的一致有界性, 此时 $0 = x_0 < x_1 < \dots < x_{J - 1} < l.$

注3.3 增长因子的计算.以Fourier展式（3.45）的通项

$$
v ^ {n} \exp (\mathrm {i} \alpha x _ {j}) \quad \left(\alpha = \frac {2 p \pi}{l}\right) \tag {3.55}
$$

代到（3.43）两端，得

$$
v ^ {n + 1} \sum_ {m \in \mathcal {N} _ {1}} a _ {m} \exp (\mathrm {i} \alpha x _ {m + j}) = v ^ {n} \sum_ {m \in \mathcal {N} _ {0}} b _ {m} \exp (\mathrm {i} \alpha x _ {m + j}).
$$

消去公因子 $\mathrm{e}^{\mathrm{i}\alpha x_j}$ ，即得

$$
v ^ {n + 1} = \left[ \sum_ {m \in \mathcal {N} _ {1}} a _ {m} \exp (\mathrm {i} \alpha x _ {m}) \right] ^ {- 1} \left[ \sum_ {m \in \mathcal {N} _ {0}} b _ {m} \exp (\mathrm {i} \alpha x _ {m}) \right] v ^ {n}.
$$

$v^n$ 前面的因子就是增长因子 $G(ph,\tau)$ （参看(3.50)).

注3.4 如果求解的是纯初值问题，则需用Fourier积分

$$
u ^ {n} (x) = \frac {1}{\sqrt {2 \pi}} \int_ {- \infty} ^ {\infty} v ^ {n} (s) \mathrm {e} ^ {\mathrm {i} x s} \mathrm {d} s \tag {3.56}
$$

代替Fourier级数(3.45).为使(3.56)有意义，应要求初值 $\varphi (x)\in L^{2}(-\infty ,\infty)$ 计算增长因子的方法如注3.3,不再详述

例3.6 考虑向前差分格式

$$
u _ {j} ^ {n + 1} = r u _ {j + 1} ^ {n} + (1 - 2 r) u _ {j} ^ {n} + r u _ {j - 1} ^ {n},
$$

以 $u_{j}^{n} = v^{n}\mathrm{e}^{\mathrm{i}\alpha jh}$ 代人，得

$$
v ^ {n + 1} \mathrm {e} ^ {\mathrm {i} \alpha j h} = \left(r \mathrm {e} ^ {\mathrm {i} \alpha (j + 1) h} + (1 - 2 r) \mathrm {e} ^ {\mathrm {i} \alpha j h} + r \mathrm {e} ^ {\mathrm {i} \alpha (j - 1) h}\right) v ^ {n}.
$$

消去 $\mathrm{e}^{\mathrm{i}\alpha jh}$ , 则知增长因子

$$
\begin{array}{l} G \left(x _ {p}, \tau\right) = (1 - 2 r) + r \left(\mathrm {e} ^ {\mathrm {i} \alpha h} + \mathrm {e} ^ {- \mathrm {i} \alpha h}\right) \\ = 1 - 2 r (1 - \cos \alpha h) \\ = 1 - 4 r \sin^ {2} \frac {\alpha h}{2}. \\ \end{array}
$$

由于 $\frac{\alpha h}{2}\left(= \frac{\pi ph}{l}\right)$ 在 $[0,\pi ]$ 中分布稠密 (随 $h\to 0$ ),为使 $G(x_{p},\tau)$ 满足von Neumann条件，必须且只需网比 $r\leqslant \frac{1}{2}$ （见3.2节例3.1)，所以向前差分格式的稳定条件是 $r\leqslant \frac{1}{2}$

注3.5Fourier方法同样可以分析差分方程组的稳定性.设差分方程组形如

$$
\sum_ {m \in \mathcal {N} _ {1}} A _ {m} U _ {j + m} ^ {n + 1} = \sum_ {m \in \mathcal {N} _ {0}} B _ {m} U _ {j + m} ^ {n}, \tag {3.57}
$$

其中 $A_{m},B_{m}$ 是 $s\times s$ 方阵，一般依赖步长 $\tau$ ，但和 $j$ 无关； $U_{j}^{n}$ 是 $s$ 维列向量,其分量为 $u_{1j}^{n},u_{2j}^{n},\dots ,u_{sj}^{n}$ .像方程式的情形一样，将 $U_{j}^{n}$ 开拓为连续变量的周期函数 $U^{n}(x) = (u_{1}^{n}(x),u_{2}^{n}(x),\dots ,u_{s}^{n}(x))^{\mathrm{T}}$ ，并将它展成Fourier级数

$$
\boldsymbol {U} ^ {n} (x) = \sum_ {p = - \infty} ^ {\infty} \boldsymbol {V} _ {p} ^ {n} \exp \left(\mathrm {i} \frac {2 p \pi}{l} x\right),
$$

将其代入 (3.57), 比较相应项的系数, 则得 $s$ 阶矩阵

$$
\boldsymbol {G} \left(x _ {p}, \tau\right) = \left[ \sum_ {m \in \mathcal {N} _ {1}} \boldsymbol {A} _ {m} \exp \left(\mathrm {i} \frac {2 m \pi}{l} x _ {p}\right) \right] ^ {- 1} \left[ \sum_ {m \in \mathcal {N} _ {0}} \boldsymbol {B} _ {m} \exp \left(\mathrm {i} \frac {2 m \pi}{l} x _ {p}\right) \right], \tag {3.58}
$$

称其为增长矩阵. 计算增长矩阵的方法跟以前一样, 以通项

$$
V ^ {n} \exp (\mathrm {i} \alpha x _ {j}) \quad \left(\alpha = \frac {2 p \pi}{l}\right) \tag {3.59}
$$

代到方程 (3.57), 消去共同因子 $\exp (\mathrm{i}\alpha x_{j})$ , 则得

$$
V ^ {n + 1} = G \left(x _ {p}, \tau\right) V ^ {n}, \tag {3.60}
$$

其中 $G(x_{p},\tau)$ 就是形如(3.58）的增长矩阵

与前面类似地有

命题3.4 差分格式(3.57)稳定的充要条件是矩阵族

$$
\left\{\boldsymbol {G} ^ {n} \left(x _ {p}, \tau\right): 0 <   \tau \leqslant \tau_ {0}, 0 <   n \tau \leqslant T, p = 0, 1, \dots , J - 1 \right\} \tag {3.61}
$$

一致有界.

命题3.5 矩阵族(3.61）一致有界的必要条件是 $G(x_{p},\tau)$ 的谱半径

$$
\rho (G) \leqslant 1 + O (\tau), \tag {3.62}
$$

即 von Neumann 条件成立.

例3.7 将Richardson格式写成等价的方程组：

$$
\left\{ \begin{array}{l} u _ {j} ^ {n + 1} = 2 r \left(u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}\right) + w _ {j} ^ {n}, \\ w _ {j} ^ {n + 1} = u _ {j} ^ {n}. \end{array} \right. \tag {3.63}
$$

以 $u_{j}^{n} = v_{1}^{n}\mathrm{e}^{\mathrm{i}\alpha x_{j}},w_{j}^{n} = v_{2}^{n}\mathrm{e}^{\mathrm{i}\alpha x_{j}}$ 代人,并消去公因子，得

$$
\left\{ \begin{array}{l} v _ {1} ^ {n + 1} = 4 r (\cos \alpha h - 1) v _ {1} ^ {n} + v _ {2} ^ {n}, \\ v _ {2} ^ {n + 1} = v _ {1} ^ {n}. \end{array} \right. \tag {3.64}
$$

显然增长矩阵

$$
\boldsymbol {G} (\alpha h) = \left[ \begin{array}{c c} - 8 r \sin^ {2} \frac {\alpha h}{2} & 1 \\ 1 & 0 \end{array} \right].
$$

由3.2节例3.4的计算结果，知 $G(\alpha h)$ 的谱半径对任意 $r > 0$ 不满足von Neumann条件，故Richardson格式恒不稳定

注3.6Fourier方法也可用于多维差分格式，求解域或为全空间（纯初值问题)，或为超长方体(周期边值条件).作为例子，考虑二维热传导方程的初边值问题

$$
\left\{ \begin{array}{l} \frac {\partial u}{\partial t} = a \left(\frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}}\right), \quad 0 <   x, y <   l (a > 0), \\ u (x, y, 0) = \varphi (x, y), \\ u (0, y, t) = u (l, y, t) = 0, \quad u (x, 0, t) = u (x, l, t) = 0. \end{array} \right. \tag {3.65}
$$

取步长 $h = \frac{l}{J},\tau = \frac{T}{N}$ ，用两族平行线 $x = x_{j} = jh,y = y_{k} = kh$ 将求解域分划成矩形网格，网点为 $(x_{j},y_{k},t_{n})(t_{n} = n\tau)$ .引进二阶差分算子

$$
\begin{array}{l} \delta_ {x} ^ {2} u _ {j k} ^ {n} = u _ {j + 1, k} ^ {n} - 2 u _ {j, k} ^ {n} + u _ {j - 1, k} ^ {n}, \\ \delta_ {y} ^ {2} u _ {j k} ^ {n} = u _ {j, k + 1} ^ {n} - 2 u _ {j, k} ^ {n} + u _ {j, k - 1} ^ {n}, \\ \end{array}
$$

作逼近 (3.65) 的向前差分格式

$$
\frac {u _ {j k} ^ {n + 1} - u _ {j k} ^ {n}}{\tau} = \frac {a}{h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j k} ^ {n} + \delta_ {y} ^ {2} u _ {j k} ^ {n}\right), \tag {3.66}
$$

向后差分格式

$$
\frac {u _ {j k} ^ {n + 1} - u _ {j k} ^ {n}}{\tau} = \frac {a}{h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j k} ^ {n + 1} + \delta_ {y} ^ {2} u _ {j k} ^ {n + 1}\right) \tag {3.67}
$$

和Crank-Nicolson格式

$$
\frac {u _ {j k} ^ {n + 1} - u _ {j k} ^ {n}}{\tau} = \frac {a}{2 h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j k} ^ {n + 1} + \delta_ {x} ^ {2} u _ {j k} ^ {n} + \delta_ {y} ^ {2} u _ {j k} ^ {n + 1} + \delta_ {y} ^ {2} u _ {j k} ^ {n}\right). \tag {3.68}
$$

它们的截断误差的阶依次为 $O\left(\tau + h^{2}\right), O\left(\tau + h^{2}\right)$ 和 $O\left(\tau^{2} + h^{2}\right)$ .

现在研究 (3.66) 的稳定性. 取通项

$$
u _ {j k} ^ {n} = v ^ {n} \exp \left(\mathrm {i} \left(\alpha x _ {j} + \beta y _ {k}\right)\right), \quad \alpha = \frac {2 \pi p}{l}, \quad \beta = \frac {2 \pi q}{l},
$$

代到 (3.66) 两端并消去公因子, 得

$$
v ^ {n + 1} = \left(1 - 4 r \sin^ {2} \frac {\alpha h}{2} - 4 r \sin^ {2} \frac {\beta h}{2}\right) v ^ {n} \quad \left(r = \frac {a \tau}{h ^ {2}}\right),
$$

从而增长因子

$$
G = G (\alpha h, \beta h) = 1 - 4 r \left(\sin^ {2} \frac {\alpha h}{2} + \sin^ {2} \frac {\beta h}{2}\right).
$$

为使 $|G| = 1 + O(\tau)$ ，必须且只需 $r \leqslant \frac{1}{4}$ . 由此可见，随着维数的增加，对网比的限制更严了.

用同样的方法可证隐格式 (3.67) 和 (3.68) 恒稳定. 逐层计算需要求解形如

$$
u _ {j k} ^ {n + 1} - c \left(\delta_ {x} ^ {2} u _ {j k} ^ {n + 1} + \delta_ {y} ^ {2} u _ {j k} ^ {n + 1}\right) = f \left(u _ {j k} ^ {n}\right) \tag {3.69}
$$

的方程 (对 (3.67), $c = r$ ; 对 (3.68), $c = \frac{r}{2}$ ), 虽然第 2 章的方法仍可以采用, 但计算量明显增加了.

# 3.3.1 习题

1. 用Fourier方法给出差分格式(3.19)的稳定条件  
2. 证明格式 (3.67) 和 (3.68) 恒稳定

# *3.4 判别差分格式稳定性的代数准则

由命题3.4知道，差分格式的稳定性归结为矩阵族(3.61)的一致有界性.当增长矩阵 $\pmb{G}$ 的阶 $s = 1$ 时， $\pmb{G}^{n}$ 一致有界的充要条件是vonNeumann条件成立，即 $|\pmb {G}| = 1 + O(\tau)$ 当 $s\geqslant 2$ 时vonNeumann条件只是稳定的必要条件，不是充分条件.但若 $\pmb {G}(x_p,\tau)$ 是正规矩阵(特别是对称矩阵)，则 $\pmb {G}(x_{p},\tau)$ 的欧氏模等于谱半径，于是vonNeumann条件也是稳定的充分条件.然而增长矩阵一般不是正规矩阵，因此有必要寻求新的充分条件.

定理3.2 设 $G(x_{p},\tau)$ 关于 $\tau$ 于 $\tau = 0$ 满足Lipschitz条件，则矩阵族(3.61）一致有界的充要条件是矩阵族

$$
\left\{G ^ {n} \left(x _ {p}, 0\right): 0 <   \tau \leqslant \tau_ {0}, \quad 0 <   n \tau \leqslant T, \quad p = 0, 1, \dots , J - 1 \right\}
$$

一致有界.

证明 只证充分性, 因为必要性的证明完全类似. 由假设, 有常数 $K$ 使 $\pmb{G}(x_p, \tau) - \pmb{G}(x_p, 0) = \tau \pmb{G}_1$ , 而 $\| \pmb{G}_1 \| \leqslant K$ . 记 $\pmb{G}_0 = \pmb{G}(x_p, 0)$ , 则 $\pmb{G} = \pmb{G}(x_p, \tau) = \pmb{G}_0 + \tau \pmb{G}_1$ .

注意

$$
\begin{array}{l} \boldsymbol {G} ^ {n} = \left(\boldsymbol {G} _ {0} + \tau \boldsymbol {G} _ {1}\right) \boldsymbol {G} ^ {n - 1} = \boldsymbol {G} _ {0} \boldsymbol {G} ^ {n - 1} + \tau \boldsymbol {G} _ {1} \boldsymbol {G} ^ {n - 1} \\ = \boldsymbol {G} _ {0} \left(\boldsymbol {G} _ {0} + \tau \boldsymbol {G} _ {1}\right) \boldsymbol {G} ^ {n - 2} + \tau \boldsymbol {G} _ {1} \boldsymbol {G} ^ {n - 1} \\ = \boldsymbol {G} _ {0} ^ {2} \boldsymbol {G} ^ {n - 2} + \tau \boldsymbol {G} _ {0} \boldsymbol {G} _ {1} \boldsymbol {G} ^ {n - 2} + \tau \boldsymbol {G} _ {1} \boldsymbol {G} ^ {n - 1} \\ \dots \\ = G _ {0} ^ {n} + \tau \sum_ {i = 0} ^ {n - 1} G _ {0} ^ {i} G _ {1} G ^ {n - i - 1}. \\ \end{array}
$$

由于 $G_0^n$ 一致有界，可设 $\| G_0^i\| \leqslant M,\| G_0^i G_1\| \leqslant M,$ ，从而得

$$
\begin{array}{l} \left\| G ^ {n} \right\| \leqslant \left\| G _ {0} ^ {n} \right\| + \tau \sum_ {i = 0} ^ {n - 1} \left\| G _ {0} ^ {i} G _ {1} \right\| \left\| G ^ {n - i - 1} \right\| \\ \leqslant M \left[ 1 + \tau \sum_ {i = 0} ^ {n - 1} \| G ^ {n - i - 1} \| \right]. \\ \end{array}
$$

由Gronwall不等式(见1.1节引理1.3)，得

$$
\left\| \boldsymbol {G} ^ {n} \right\| \leqslant M (1 + M \tau) ^ {n - 1}.
$$

又 $0 < n \leqslant \frac{T}{\tau}$ , 故不等式右端一致有界.

例3.8 考虑逼近带低阶项的抛物方程：

$$
\frac {\partial u}{\partial t} = \frac {\partial^ {2} u}{\partial x ^ {2}} + b u.
$$

逼近它的向前差分格式为

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} + b u _ {j} ^ {n}.
$$

用Fourier方法即知增长因子 $G = 1 - 4r\sin^2{\frac{\alpha h}{2}} + b\tau.$ 由定理3.2, $G^n$ 一致有界等价于 $G_0^n\left(G_0 = 1 - 4r\sin^2{\frac{\alpha h}{2}}\right)$ 一致有界.由例3.6, $G_0^n$ 一致有界的充要条件是 $r\leqslant \frac{1}{2}$ 这说明低阶项 $bu_{j}^{n}$ 不影响差分格式的稳定性

今设 $G(x_{p},\tau) = G(x_{p})$ 与 $\tau$ 无关.考虑矩阵族

$$
\left\{\boldsymbol {G} ^ {n} \left(x _ {p}\right): x _ {0} = 0 <   x _ {1} <   \dots <   x _ {J} = l, \quad 0 <   \tau \leqslant \tau_ {0}, \quad 0 <   n \tau \leqslant T \right\} \tag {3.70}
$$

的一致有界性.

命题3.6 矩阵族(3.70)一致有界的充要条件是矩阵族

$$
\left\{\boldsymbol {G} ^ {n} (x): 0 \leqslant x \leqslant l, \quad n = 1, 2, \dots \right\} \tag {3.71}
$$

一致有界.

证明 充分性显然, 只证必要性. 假定网格按 2 等分, 4 等分, $\cdots, 2^{m}$ 等分加密, 则二等分点 $x_{j}$ 一旦是 $[0, l]$ 的网点便永远是网点. 由假设,

$$
\| G ^ {n} (x _ {j}) \| \leqslant M, \quad 0 <   n \tau \leqslant T,
$$

$M$ 是与分划无关的常数. 令 $\tau \to 0$ (从而 $h \to 0$ ), 则

$$
\| G ^ {n} \left(x _ {j}\right) \| \leqslant M, \quad n = 1, 2, \dots ,
$$

而二等分点 $\{x_{j}\}$ 于 $[0,l]$ 稠密, $G(x)$ 是连续函数, 故

$$
\| G ^ {n} (x) \| \leqslant M, \quad 0 \leqslant x \leqslant l, \quad n = 1, 2, \dots .
$$

定理3.3（一致对角化） 若对任一 $G(x_{p},\tau)$ ，有矩阵 $\pmb{H}$ 使

$$
\boldsymbol {H} ^ {- 1} \boldsymbol {G} \boldsymbol {H} = \boldsymbol {\Lambda} = \left[ \begin{array}{c c c} \lambda_ {1} & & \\ & \ddots & \\ & & \lambda_ {s} \end{array} \right],
$$

且 $\pmb{H}$ 和 $H^{-1}$ 关于 $p$ 和充分小的 $\tau >0$ 一致有界，则von Neumann条件也是稳定的充分条件.

证明 因 $G^{n} = H\Lambda^{n}H^{-1}$ ，故结论显然

若差分方程组 (3.57) 是二阶方程组, 则增长矩阵 $G(x, \tau)$ 也是二阶, 此时有一个便于检验的稳定性条件.

定理3.4 设 $G(x, \tau)$ 是二阶矩阵， $g_{ij}$ 是 $G$ 的第 $i$ 行第 $j$ 列元素， $\lambda_1$ 和 $\lambda_2$ 是 $G$ 的特征值，若下列条件成立：

$(\alpha)$ $|\lambda_{i}(x,\tau)|\leqslant 1 + M\tau ,\quad i = 1,2,$

$(\beta)$ $||G(x,\tau) - \frac{1}{2} (g_{11}(x,\tau) + g_{22}(x,\tau))I||$

$$
\leqslant M (\tau + | 1 - | \lambda_ {1} (x, \tau) | | + | \lambda_ {1} (x, \tau) - \lambda_ {2} (x, \tau) |),
$$

则矩阵族 (3.61) 一致有界, 其中 $\pmb{I}$ 是二阶单位矩阵.

证明 注意条件 $(\beta)$ 实际上是

$$
\left\| \boldsymbol {G} - \frac {1}{2} \left(\lambda_ {1} + \lambda_ {2}\right) \boldsymbol {I} \right\| \leqslant M \left(\tau + | 1 - | \lambda_ {1} | | + | \lambda_ {1} - \lambda_ {2} |\right). \tag {3.72}
$$

利用 $G$ 的Jordan标准形，我们有

$$
\boldsymbol {G} ^ {n} - \lambda^ {n} \boldsymbol {I} = n \lambda^ {n - 1} (\boldsymbol {G} - \lambda \boldsymbol {I}), \quad \lambda_ {1} = \lambda_ {2} = \lambda , \tag {3.73}
$$

$$
\boldsymbol {G} ^ {n} - \frac {1}{2} \left(\lambda_ {1} ^ {n} + \lambda_ {2} ^ {n}\right) \boldsymbol {I} = \frac {\lambda_ {1} ^ {n} - \lambda_ {2} ^ {n}}{\lambda_ {1} - \lambda_ {2}} \left(\boldsymbol {G} - \frac {1}{2} \left(\lambda_ {1} + \lambda_ {2}\right) \boldsymbol {I}\right), \quad \lambda_ {1} \neq \lambda_ {2}. \tag {3.74}
$$

形式上也可将(3.73)写成(3.74).设 $\lambda_{1} = r_{1}\mathrm{e}^{\mathrm{i}\theta_{1}},\lambda_{2} = r_{2}\mathrm{e}^{\mathrm{i}\theta_{2}}$ ，不妨设 $r_1\geqslant r_2$

由条件 $(\alpha), \lambda_1^n$ 和 $\lambda_2^n$ 一致有界. 由 (3.74), 可将 $G^n$ 一致有界归结为 (3.74) 右端一致有界. 因 $\frac{|\lambda_1^n - \lambda_2^n|}{|\lambda_1 - \lambda_2|} \leqslant nr_1^{n-1}, n\tau \leqslant T,$ 由 (3.72) 可知只需证明 $nr_1^{n-1}|1-r_1|$ 一致有界. 当 $r_1 < 1$ 时, $nr_1^{n-1}(1-r_1) \leqslant (1+r_1+\cdots+r_1^{n-1})(1-r_1)=1-r_1^n \leqslant 1$ ; 当 $r_1 \geqslant 1$ 时由 $(\alpha)$ 知 $|1-r_1| \leqslant M\tau$ , 从而 $nr_1^{n-1}|1-r_1| \leqslant Mn\tau \leqslant M_T$ . □

注3.7 可以证明， $(\alpha),(\beta)$ 也是矩阵族(3.61）一致有界的必要条件（参见[4]).

推论3.3 特别地，若 $G(x,\tau)$ 与 $\tau$ 无关，则知二阶矩阵族(3.71）一致有界的充要条件是

$$
\begin{array}{l} \left(\alpha\right) ^ {\prime} \quad \left| \lambda_ {1} (x) \right| \leqslant 1, \quad i = 1, 2, \quad 0 \leqslant x \leqslant l, \\ \left(\beta\right) ^ {\prime} \quad \| \boldsymbol {G} (x) - \frac {1}{2} \left(g _ {1 1} (x) + g _ {2 2} (x)\right) \boldsymbol {I} \| \\ \leqslant M \left(| 1 - | \lambda_ {1} (x) | | + | \lambda_ {1} (x) - \lambda_ {2} (x) |\right), \quad 0 \leqslant x \leqslant l. \\ \end{array}
$$

注3.8 定理3.4是一个很有用的稳定性判别法，应用时注意到以下两点是方便的：

1. 实系数二次方程 $\lambda^2 - b\lambda - c = 0$ 的根按模不大于 1 的充要条件是

$$
| b | \leqslant 1 - c \leqslant 2. \tag {3.75}
$$

这在检验 von Neumann 条件时有用.

2. 检验条件 $(\beta)$ 时要计算二阶矩阵的范数，通常用 Frobenius 范数 $\| \cdot \|_{\mathrm{F}}$ (F 范数). $m$ 阶矩阵 $\mathbf{A} = (a_{ij})$ 的 F 范数是

$$
\| \boldsymbol {A} \| _ {\mathrm {F}} = \left(\sum_ {i, j = 1} ^ {m} | a _ {i j} | ^ {2}\right) ^ {\frac {1}{2}}. \tag {3.76}
$$

例3.9 考虑逼近热传导方程的Dufort-Frankel格式

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{2 \tau} = a \frac {u _ {j + 1} ^ {n} - u _ {j} ^ {n + 1} - u _ {j} ^ {n - 1} + u _ {j - 1} ^ {n}}{h ^ {2}}. \tag {3.77}
$$

引进新变量 $v_{j}^{n + 1} = u_{j}^{n}$ ，将它化为一阶方程组：

$$
\left\{ \begin{array}{l} u _ {j} ^ {n + 1} = \frac {2 r}{1 + 2 r} \left(u _ {j + 1} ^ {n} + u _ {j - 1} ^ {n}\right) + \frac {1 - 2 r}{1 + 2 r} v _ {j} ^ {n}, \\ v _ {j} ^ {n + 1} = u _ {j} ^ {n}, \quad r = \frac {\tau}{h ^ {2}} \end{array} \right.
$$

或

$$
\binom {u _ {j} ^ {n + 1}} {v _ {j} ^ {n + 1}} = \left[ \begin{array}{c c} \frac {2 r}{1 + 2 r} (T _ {1} + T _ {- 1}) & \frac {1 - 2 r}{1 + 2 r} \\ 1 & 0 \end{array} \right] \binom {u _ {j} ^ {n}} {v _ {j} ^ {n}},
$$

其中 $T_{\pm 1}u_{j} = u_{j\pm 1}$ 为移位算子.用Fourier方法可知增长矩阵

$$
\boldsymbol {G} (\alpha h) = \left[ \begin{array}{c c} \frac {2 r}{1 + 2 r} \left(\mathrm {e} ^ {\mathrm {i} \alpha h} + \mathrm {e} ^ {- \mathrm {i} \alpha h}\right) & \frac {1 - 2 r}{1 + 2 r} \\ 1 & 0 \end{array} \right] = \left[ \begin{array}{c c} \frac {4 r}{1 + 2 r} \cos \alpha h & \frac {1 - 2 r}{1 + 2 r} \\ 1 & 0 \end{array} \right].
$$

其特征方程为

$$
\lambda^ {2} - \frac {4 r \cos \theta}{1 + 2 r} \lambda - \frac {1 - 2 r}{1 + 2 r} = 0, \quad \theta = \alpha h = \frac {2 \pi p h}{l} \in [ 0, 2 \pi ], \tag {3.78}
$$

它的系数显然满足条件(3.75)，故其特征值按模 $\leqslant 1$ ，从而条件 $(\alpha)^{\prime}$ 成立

其次，(3.78) 的二根为

$$
\lambda_ {1, 2} = \frac {2 r \cos \theta \pm \sqrt {1 - 4 r ^ {2} \sin^ {2} \theta}}{1 + 2 r}. \tag {3.79}
$$

令 $\varLambda(\theta)=|1-|\lambda_1||+|\lambda_1-\lambda_2|$ , 其中

$$
| \lambda_ {1} - \lambda_ {2} | = \frac {2 \left| \sqrt {1 - 4 r ^ {2} \sin^ {2} \theta} \right|}{1 + 2 r}.
$$

当 $\lambda_1, \lambda_2$ 是实根时， $|\lambda_1| = \frac{2r|\cos\theta|}{1 + 2r} + \frac{\left|\sqrt{1 - 4r^2\sin^2\theta}\right|}{1 + 2r}$ ，所以

$$
\begin{array}{l} \Lambda (\theta) = 1 - \frac {2 r | \cos \theta |}{1 + 2 r} + \frac {\left| \sqrt {1 - 4 r ^ {2} \sin^ {2} \theta} \right|}{1 + 2 r} \geqslant 1 - \frac {2 r}{1 + 2 r} \\ = \frac {1}{1 + 2 r} > 0. \\ \end{array}
$$

当 $\lambda_1, \lambda_2$ 是复根时， $\lambda_2 = \bar{\lambda}_1, |\lambda_1| = \sqrt{|\lambda_1\lambda_2|} = \sqrt{\left|\frac{1 - 2r}{1 + 2r}\right|}$ ，从而

$$
\Lambda (\theta) \geqslant 1 - | \lambda_ {1} | = 1 - \sqrt {\frac {| 1 - 2 r |}{| 1 + 2 r |}} > 0.
$$

可见函数 $\varLambda(\theta)$ 对任意 $r > 0$ 于 $[0,2\pi ]$ 上有正的下界 $m > 0$

另一方面，

$$
\boldsymbol {G} (\theta) - \frac {1}{2} \left(g _ {1 1} + g _ {2 2}\right) \boldsymbol {I} = \left[ \begin{array}{c c} \frac {2 r}{1 + 2 r} \cos \theta & \frac {1 - 2 r}{1 + 2 r} \\ 1 & - \frac {2 r}{1 + 2 r} \cos \theta \end{array} \right]
$$

其F范数显然有上界 $K > 0$ ，故条件 $(\beta)'$ 成立.由定理3.4的推论，可知(3.77） $\forall r > 0$ 都稳定.

判断差分格式稳定性最完整的代数准则已由Kreiss得到[28],但这些准则难以检验,所以有不少工作研究便于应用的充分条件(参见[14,28]).

# 3.4.1 习题

1. 证明实系数二次方程 $\lambda^2 - b\lambda - c = 0$ 的根按模小于或等于1的充要条件是

$$
| b | \leqslant 1 - c \leqslant 2.
$$

2. 证明差分格式 (3.21) 恒稳定  
3. 证明差分格式

$$
\left\{ \begin{array}{l} u _ {j} ^ {n + 1} - u _ {j} ^ {n} = a r \left(u _ {j + 1} ^ {n} - u _ {j} ^ {n} - u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}\right), \\ u _ {j} ^ {n + 2} - u _ {j} ^ {n + 1} = a r \left(u _ {j + 1} ^ {n + 2} - u _ {j} ^ {n + 2} - u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}\right) \quad (a > 0) \end{array} \right.
$$

(Saul'yev, 1957) 恒稳定.

4. 考查如下隐-显格式：

$$
\left\{ \begin{array}{l} \frac {u _ {j} ^ {2 m - 1} - u _ {j} ^ {2 m - 2}}{\tau} = a \frac {u _ {j + 1} ^ {2 m - 2} - 2 u _ {j} ^ {2 m - 2} + u _ {j - 1} ^ {2 m - 2}}{h ^ {2}}, \\ \frac {u _ {j} ^ {2 m} - u _ {j} ^ {2 m - 1}}{\tau} = a \frac {u _ {j + 1} ^ {2 m} - 2 u _ {j} ^ {2 m} + u _ {j - 1} ^ {2 m}}{h ^ {2}}, \\ u _ {j} ^ {0} = \phi_ {j}, \quad j = 1, 2, \dots , J - 1, \quad m = 1, 2, \dots , M, \\ u _ {0} = u _ {J} = 0 \end{array} \right. (a > 0)
$$

其计算量大约是向后差分格式的一半. 试证明它是恒稳定的.

# 第4章

# 双曲型方程的有限差分法

# 4.1 波动方程的差分逼近

![](images/c0229b1cf12dc9eb3766a3fe1b01a30ddc06366977a9976525413b5b1e7472b0.jpg)

# 人物简介

# 4.1.1 波动方程及其特征

二阶线性双曲型偏微分方程的最简单模型是波动方程：

$$
\frac {\partial^ {2} u}{\partial t ^ {2}} = a ^ {2} \frac {\partial^ {2} u}{\partial x ^ {2}}, \tag {4.1}
$$

其中 $a > 0$ 是常数. 根据二阶偏微分方程理论, 与 (4.1) 相应的特征方程为

$$
\mathrm {d} x ^ {2} - a ^ {2} \mathrm {d} t ^ {2} = 0
$$

或

$$
1 - a ^ {2} \left(\frac {\mathrm {d} t}{\mathrm {d} x}\right) ^ {2} = 0. \tag {4.2}
$$

由此定出两个方向：

$$
\frac {\mathrm {d} t}{\mathrm {d} x} = \pm \frac {1}{a}, \tag {4.3}
$$

称为特征方向. 解常微分方程 (4.3), 得到两族直线:

$$
x - a t = c _ {1}, \quad x + a t = c _ {2}, \tag {4.4}
$$

称为特征.

在研究波动方程的各种定解问题时，特征起着重要作用。例如，我们用 $u$ 沿特征的偏导数表示它沿 $x, t$ 的偏导数，则

$$
\frac {\partial^ {2} u}{\partial t ^ {2}} = a ^ {2} \left(\frac {\partial^ {2} u}{\partial c _ {1} ^ {2}} - 2 \frac {\partial^ {2} u}{\partial c _ {1} \partial c _ {2}} + \frac {\partial^ {2} u}{\partial c _ {2} ^ {2}}\right),
$$

$$
\frac {\partial^ {2} u}{\partial x ^ {2}} = \frac {\partial^ {2} u}{\partial c _ {1} ^ {2}} + 2 \frac {\partial^ {2} u}{\partial c _ {1} \partial c _ {2}} + \frac {\partial^ {2} u}{\partial c _ {2} ^ {2}},
$$

于是方程（4.1）化为

$$
\frac {\partial^ {2} u}{\partial c _ {1} \partial c _ {2}} = 0,
$$

从而其通解为

$$
u = f _ {1} \left(c _ {1}\right) + f _ {2} \left(c _ {2}\right) = f _ {1} (x - a t) + f _ {2} (x + a t).
$$

如果 $u$ 在 $x$ 轴的初值为

$$
u (x, 0) = \varphi_ {0} (x), u _ {t} (x, 0) = \varphi_ {1} (x), - \infty <   x <   \infty , \tag {4.5}
$$

则可以定出 $f_{1}, f_{2}$ , 从而得

$$
u (x, t) = \frac {1}{2} \left[ \varphi_ {0} (x + a t) + \varphi_ {0} (x - a t) \right] + \frac {1}{2 a} \int_ {x - a t} ^ {x + a t} \varphi_ {1} (\xi) \mathrm {d} \xi . \tag {4.6}
$$

这是熟知的 d'Alembert 公式

公式（4.6）告诉我们， $u$ 在点 $(x_0,t_0)(t_0 > 0)$ 的值仅依赖于初值函数 $\varphi_0(x),\varphi_1(x)$ 在区间 $[x_0 - at_0,x_0 + at_0]$ 上的值，与区间外的初值无关，故称 $[x_0 - at_0,x + at_0]$ 为点 $(x_0,t_0)$ 的依存域.其实，区间 $[x_0 - at_0,x_0 + at_0]$ 上的初值不只确定了 $u(x_0,t_0)$ ，而且确定了 $u$ 在以 $(x_0 - at_0,0),(x_0 + at_0,0),(x_0,t_0)$ 为顶点的三角形域内的值，故称此三角形域为区间 $[x_0 - at_0,x_0 + at_0]$ 的决定域

显然, 为了得到点 $(x_0, t_0)$ 的依存域, 只需通过 $(x_0, t_0)$ 作两条特征, 它们与 $x$ 轴截出的闭区间即是. 为了得到区间 $[x_0 - at_0, x_0 + at_0]$ 的决定域, 过 $(x_0 - at_0, 0)$ 作第一特征 (斜率为正), 过 $(x_0 + at_0, 0)$ 作第二特征 (斜率为负), 它们交出的三角形域即为决定域 (参看图 4.1(a)). 从公式 (4.6) 还看出, 对 $x$ 轴上任一点 $(x_0, 0)$ , 依存域包含 $(x_0, 0)$ 的一切 $(x, t)$ 的集合恰好是以 $(x_0, 0)$ 为顶点, 过 $(x_0, 0)$ 的特征 $x - at = x_0$ 和 $x + at = x_0 (t > 0)$ 为边的角形域, 称之为 $(x_0, 0)$ 的影响域 (参看图 4.1(b)).

![](images/40269e032ed2cd6ac19ebe40a619535ea7083322804f1230c720f98efff3a347.jpg)  
(a) 决定域

![](images/7435ccb8d175d7b2c5c98f1837b126502fc9c1c2ae1211fbd8eaab57a8b551b7.jpg)  
(b) 影响域  
图4.1 决定域和影响域

# 4.1.2 显格式

现在构造 (4.1) 的差分逼近. 取空间步长 $h$ 和时间步长 $\tau$ , 用两族平行直线

$$
x = x _ {j} = j h, \quad j = 0, \pm 1, \pm 2, \dots ,
$$

$$
t = t _ {n} = n \tau , \quad n = 0, 1, 2 \dots
$$

作矩形网格. 于网点 $(x_{j}, t_{n})$ 用 Taylor 展开式, 得

$$
\frac {u \left(x _ {j + 1} , t _ {n}\right) - 2 u \left(x _ {j} , t _ {n}\right) + u \left(x _ {j - 1} , t _ {n}\right)}{h ^ {2}} = u _ {x x} ^ {\prime \prime} \left(x _ {j}, t _ {n}\right) + \frac {h ^ {2}}{1 2} \frac {\partial^ {4}}{\partial x ^ {4}} u \left(x _ {j}, t _ {n}\right) + O \left(h ^ {4}\right),
$$

$$
\frac {u \left(x _ {j} , t _ {n + 1}\right) - 2 u \left(x _ {j} , t _ {n}\right) + \left(x _ {j} , t _ {n - 1}\right)}{\tau^ {2}} = u _ {t t} ^ {\prime \prime} \left(x _ {j}, t _ {n}\right) + \frac {\tau^ {2}}{1 2} \frac {\partial^ {4}}{\partial t ^ {4}} u \left(x _ {j}, t _ {n}\right) + O \left(\tau^ {4}\right),
$$

将 $u_{xx}^{\prime \prime}(x_j,t_n),u_{tt}^{\prime \prime}(x_j,t_n)$ 代入波动方程，得

$$
\begin{array}{l} \frac {u \left(x _ {j} , t _ {n + 1}\right) - 2 u \left(x _ {j} , t _ {n}\right) + \left(x _ {j} , t _ {n - 1}\right)}{\tau^ {2}} \\ = a ^ {2} \frac {u (x _ {j + 1} , t _ {n}) - 2 u (x _ {j} , t _ {n}) + u (x _ {j - 1} , t _ {n})}{h ^ {2}} + R _ {j} ^ {n} (u), \\ \end{array}
$$

其中

$$
R _ {j} ^ {n} (u) = \frac {h ^ {2}}{1 2 a ^ {2}} \left(r ^ {2} \frac {\partial^ {4}}{\partial t ^ {4}} u \left(x _ {j}, t _ {n}\right) - a ^ {4} \frac {\partial^ {4}}{\partial x ^ {4}} u \left(x _ {j}, t _ {n}\right)\right) + O \left(\tau^ {4} + h ^ {4}\right) \tag {4.7}
$$

是截断误差, $r = \frac{ah}{\tau}$ 是网比. 舍去 $R_{j}^{n}(u)$ , 则得差分方程:

$$
\frac {u _ {j} ^ {n + 1} - 2 u _ {j} ^ {n} + u _ {j} ^ {n - 1}}{\tau^ {2}} = a ^ {2} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}, \quad j = 0, \pm 1, \pm 2, \dots , \quad n = 1, 2, \dots , \tag {4.8}
$$

这里 $u_{j}^{n}$ 表示 $\pmb{u}$ 于 $(x_{j},t_{n})$ 的近似值.初值条件用下列差分方程代替：

$$
u _ {j} ^ {0} = \varphi_ {0} \left(x _ {j}\right), \tag {4.9}
$$

$$
\frac {u _ {j} ^ {1} - u _ {j} ^ {0}}{\tau} = \varphi_ {1} \left(x _ {j}\right). \tag {4.10}
$$

显然，(4.8)的截断误差的阶是 $O\left(\tau^2 + h^2\right)$ ，而(4.10)的截断误差的阶仅为 $O(\tau)$ . 为了提高精度，也可以用中心差商代替 $u_t'$ ，得

$$
\frac {u _ {j} ^ {1} - u _ {j} ^ {- 1}}{2 \tau} = \varphi_ {1} \left(x _ {j}\right). \tag {4.11}
$$

于 (4.8) 令 $n = 0$ ，又得

$$
\frac {u _ {j} ^ {1} - 2 u _ {j} ^ {0} + u _ {j} ^ {- 1}}{\tau^ {2}} = a ^ {2} \frac {u _ {j + 1} ^ {0} - 2 u _ {j} ^ {0} + u _ {j - 1} ^ {0}}{h ^ {2}}.
$$

消去 $u_{j}^{-1}$ ，则

$$
u _ {j} ^ {1} = \frac {r ^ {2}}{2} \left[ \varphi_ {0} \left(x _ {j - 1}\right) + \varphi_ {0} \left(x _ {j + 1}\right) \right] + (1 - r ^ {2}) \varphi_ {0} \left(x _ {j}\right) + \tau \varphi_ {1} \left(x _ {j}\right). \tag {4.12}
$$

利用 (4.9) (4.10) (或 (4.12)) 可算出初始层 $(n = 0)$ 及第一层 $(n = 1)$ 各网格节点上的值. 然后利用 (4.8) 或

$$
u _ {j} ^ {n + 1} = r ^ {2} \left(u _ {j - 1} ^ {n} + u _ {j + 1} ^ {n}\right) + 2 \left(1 - r ^ {2}\right) u _ {j} ^ {n} - u _ {j} ^ {n - 1}, \tag {4.13}
$$

就可逐层算出任意网点的值

公式 (4.8) 是显式的三层差分格式, 节点分布如图 4.2. 上述格式也可用以解混合问题

$$
\left\{ \begin{array}{l l} \frac {\partial^ {2} u}{\partial t ^ {2}} = a ^ {2} \frac {\partial^ {2} u}{\partial x ^ {2}}, & 0 <   x <   l, \quad 0 <   t \leqslant T, \\ u (x, 0) = \varphi_ {0} (x), & u _ {t} ^ {\prime} (x, 0) = \varphi_ {1} (x), \\ u (0, t) = \alpha (t), & u (l, t) = \beta (t). \end{array} \right. \tag {4.14}
$$

这时取 $h = \frac{l}{J},\tau = \frac{T}{N}$ 除(4.8)—(4.10)外，再补充边值条件

$$
u _ {0} ^ {n} = \alpha (n \tau), \quad u _ {J} ^ {n} = \beta (n \tau). \tag {4.15}
$$

![](images/53daf5e96714bb21c9b31144a99149ccb077cea9fc9f191b7e1652b24ebeba40.jpg)  
图4.2 节点分布

# 4.1.3 稳定性分析

像抛物型方程一样，造出差分格式后，要检验它是否稳定，在什么条件下稳定，这是差分方程理论的基本问题.

为了引用第3章判别稳定性的方法，我们把波动方程（4.1）化成一阶微分方程组，相应地把三层差分格式（4.8）化成两层差分格式。一种简单的做法是引进变量 $v = \frac{\partial u}{\partial t}$ ，于是（4.1）化为

$$
\frac {\partial u}{\partial t} = v, \quad \frac {\partial v}{\partial t} = a ^ {2} \frac {\partial^ {2} u}{\partial x ^ {2}}.
$$

但对于构造差分逼近，更常用的方法是再引进变量 $w = a\frac{\partial u}{\partial x}$ ，将（4.1）化为

$$
\frac {\partial v}{\partial t} = a \frac {\partial w}{\partial x}, \quad \frac {\partial w}{\partial t} = a \frac {\partial v}{\partial x}. \tag {4.16}
$$

若令 $\pmb {U} = (v,w)^{\mathrm{T}}$

$$
\boldsymbol {A} = \left[ \begin{array}{c c} 0 & a \\ a & 0 \end{array} \right],
$$

则 (4.16) 可写成为

$$
\frac {\partial \boldsymbol {U}}{\partial t} - \boldsymbol {A} \frac {\partial \boldsymbol {U}}{\partial x} = \mathbf {0}. \tag {4.17}
$$

相应地，将(4.8)写成等价的双层差分格式：

$$
\left\{ \begin{array}{l} \frac {v _ {j} ^ {n + 1} - v _ {j} ^ {n}}{\tau} = a \frac {w _ {j + \frac {1}{2}} ^ {n} - w _ {j - \frac {1}{2}} ^ {n}}{h}, \\ \frac {w _ {j - \frac {1}{2}} ^ {n + 1} - w _ {j - \frac {1}{2}} ^ {n}}{\tau} = a \frac {v _ {j} ^ {n + 1} - v _ {j - 1} ^ {n + 1}}{h}, \end{array} \right. \tag {4.18}
$$

其中

$$
v _ {j} ^ {n} = \frac {u _ {j} ^ {n} - u _ {j} ^ {n - 1}}{\tau}, \quad w _ {j - \frac {1}{2}} ^ {n + 1} = a \frac {u _ {j} ^ {n + 1} - u _ {j - 1} ^ {n + 1}}{h}.
$$

现在用Fourier方法分析（4.18）的稳定性.为此，考虑具周期边值条件的混合问题.按照第3章第3节，以

$$
v _ {j} ^ {n} = V _ {1} ^ {n} \exp (\mathrm {i} \alpha x _ {j}),
$$

$$
w _ {j} ^ {n} = V _ {2} ^ {n} \exp (\mathrm {i} \alpha x _ {j})
$$

$\left(\alpha = \frac{2p\pi}{l},\mathrm{i} = \sqrt{-1}\right)$ 代到(4.17)，消去共同因子 $\exp (\mathrm{i}\alpha x_j)$ 和 $\exp \left(\mathrm{i}\alpha x_{j - \frac{1}{2}}\right)$ ，得

$$
\begin{array}{l} V _ {1} ^ {n + 1} - 2 \mathrm {i} r \left(\sin \frac {\pi p h}{l}\right) V _ {2} ^ {n} = V _ {1} ^ {n}, \\ - 2 \mathrm {i} r \left(\sin \frac {\pi p h}{l}\right) V _ {1} ^ {n + 1} + V _ {2} ^ {n + 1} = V _ {2} ^ {n}, \\ \end{array}
$$

或

$$
\left( \begin{array}{c} V _ {1} ^ {n + 1} \\ V _ {2} ^ {n + 1} \end{array} \right) = G \left(\frac {\pi p h}{l}\right) \left( \begin{array}{c} V _ {1} ^ {n} \\ V _ {2} ^ {n} \end{array} \right),
$$

其中

$$
\boldsymbol {G} \left(\frac {\pi p h}{l}\right) = \left[ \begin{array}{c c} 1 & \mathrm {i} c \\ \mathrm {i} c & 1 - c ^ {2} \end{array} \right] \left(c = 2 r \sin \frac {\pi p h}{l}\right) \tag {4.19}
$$

为增长矩阵, $r = \frac{a\tau}{h}$ 为网比.由命题3.6，差分格式(4.18)稳定的充要条件是矩阵族

$$
\left\{\boldsymbol {G} ^ {n} (\theta) \right\} (0 \leqslant \theta \leqslant \pi , n = 1, 2 \dots)
$$

一致有界.

注意 $G(\theta)$ 的特征方程为

$$
\lambda^ {2} - (2 - c ^ {2}) \lambda + 1 = 0, \tag {4.20}
$$

它的根按模 $\leqslant 1$ 的充要条件是（见3.4.1小节习题1）

$$
\left| 2 - c ^ {2} \right| \leqslant 2,
$$

即 $r \leqslant 1$ . 这是差分格式稳定的必要条件——von Neumann 条件.

方程（4.20）的二根为

$$
\lambda_ {1, 2} = \frac {1}{2} \left[ (2 - c ^ {2}) \pm \mathrm {i} | c | \sqrt {4 - c ^ {2}} \right],
$$

故

$$
\left| \lambda_ {1} - \lambda_ {2} \right| = | c | \sqrt {4 - c ^ {2}}.
$$

又因 $c^2 = 4r^2\sin^2\theta \leqslant 4$ 故 $\lambda_2 = \bar{\lambda}_1,|\lambda_1| = |\lambda_2| = 1,$ 从而

$$
1 - | \lambda_ {1} | = 0.
$$

另一方面，

$$
\boldsymbol {G} - \frac {1}{2} \left(\lambda_ {1} + \lambda_ {2}\right) \boldsymbol {I} = \left[ \begin{array}{c c} \frac {1}{2} c ^ {2} & \mathrm {i} c \\ \mathrm {i} c & - \frac {1}{2} c ^ {2} \end{array} \right],
$$

其F-模为

$$
\left\| \boldsymbol {G} - \frac {1}{2} \left(\lambda_ {1} + \lambda_ {2}\right) \boldsymbol {I} \right\| _ {\mathrm {F}} = | c | \left(2 + \frac {1}{2} c ^ {2}\right) ^ {\frac {1}{2}}.
$$

为使定理3.4推论中的条件 $(\beta)$ 成立，只需存在 $M > 0$ 使

$$
\left| c \right| \left(2 + \frac {1}{2} c ^ {2}\right) ^ {\frac {1}{2}} \leqslant M \left| \lambda_ {1} - \lambda_ {2} \right| = M \left| c \right| \sqrt {4 - c ^ {2}}. \tag {4.21}
$$

上式当 $r < 1$ 时显然成立. 若 $r = 1$ , 则当 $\theta = \frac{\pi}{2}$ 时, $c^2 = 4$ , 上式右端为0, 不等式(4.21)不成立, 故应要求 $r < 1$ .

总之，差分格式(4.18)稳定的充要条件是网比

$$
r = \frac {a \tau}{h} <   1. \tag {4.22}
$$

注4.1 其实 $r = 1$ 时, $\pmb{G}(\theta)$ 于 $\theta = \frac{\pi}{2}$ 有重根 $\lambda_1 = \lambda_2 = -1$ , 且初等因子的次数等于2, 因此有相似变换 $\pmb{S}$ , 使

$$
\boldsymbol {G} \left(\frac {\pi}{2}\right) = \boldsymbol {S} \left[ \begin{array}{c c} - 1 & 1 \\ 0 & - 1 \end{array} \right] \boldsymbol {S} ^ {- 1},
$$

从而

$$
\boldsymbol {G} ^ {n} \left(\frac {\pi}{2}\right) = \boldsymbol {S} \left[ \begin{array}{c c} (- 1) ^ {n} & (- 1) ^ {n - 1} n \\ 0 & (- 1) ^ {n} \end{array} \right] \boldsymbol {S} ^ {- 1}, \quad n = 1, 2, \dots ,
$$

这表明 $G^{n}\left(\frac{\pi}{2}\right)$ 无界，所以当 $r = 1$ 时格式(4.18)不稳定，但它关于 $n$ 是线性增长，所以也称线性不稳定。如果方程的解充分光滑，则差分格式截断误差的阶达到2。按定理3.1的证法仍可证明差分解收敛（参看[31])。

稳定性条件 (4.22) 有一直观几何解释. 从方程 (4.13) 看出, $u_{j}^{n}$ 依赖前两层值 $u_{j-1}^{n-1}, u_{j}^{n-1}, u_{j+1}^{n-1}, u_{j}^{n-2}$ , 这四个值又依赖 $u_{j-2}^{n-2}, u_{j-1}^{n-2}, u_{j}^{n-2}, u_{j+1}^{n-2}, u_{j+2}^{n-2}$ 和 $u_{j-1}^{n-3}, u_{j}^{n-3}, u_{j+1}^{n-3}$ . 以此类推, 可知 $u_{j}^{n}$ 最终依赖初始层 $n = 0$ 上的下列值:

$$
u _ {j - n} ^ {0}, u _ {j - n + 1} ^ {0}, \dots , u _ {j} ^ {0}, \dots , u _ {j + n - 1} ^ {0}, u _ {j + n} ^ {0}.
$$

因此称 $x$ 轴上含于区间 $[x_{j - n}, x_{j + n}]$ 的网点为差分解 $u_{j}^{n}$ 的依存域，它是 $x$ 轴上被过 $(x_{j}, t_{n})$ 的两条直线

$$
x - x _ {j} = \pm \frac {h}{\tau} (t - t _ {n})
$$

截下的区间所覆盖的网域. 注意过 $(x_{j}, t_{n})$ 的两条特征线为 $x - x_{j} = \pm a (t - t_{n})$ , 差分方程稳定性的必要条件为 $r \leqslant 1$ 或

$$
\frac {\tau}{h} \leqslant a ^ {- 1}.
$$

![](images/1e17bccac91dd9ed45ccba5579f4fdd93d85ce270eb1c083eb6b4bf2e694224b.jpg)

人物简介

可见差分方程稳定的必要条件是差分解的依存域必须包含微分方程解的依存域，否则差分方程不稳定.

现在利用依存域的概念证明：当 $r > 1$ 时差分解不收敛。如图4.3， $r > 1$ ，微分方程解的依存域 $[P', Q']$ 大于差分解的依存域 $[P, Q]$ 。固定 $(x_j, t_n)$ ，让网格步长变小，但网比 $r$ 保持不变，则依存域 $[P', Q']$ ， $[P, Q]$ 不变。显然，若改变区间 $(P', P)$ 和 $(Q, Q')$ 上的初值，但 $[P, Q]$ 上的初值不变，则 $u(x_j, t_n)$ 可取不同值，而 $u_j^n$ 当 $h \to 0, \tau \to 0$ 时（ $r$ 固定不变）是一串确定的数列，它不可能收敛到不同的 $u(x_j, t_n)$ 。总之我们知道，当 $r < 1$ 时，差分方程稳定，因而差分解收敛（参看3.2节）。Courant等曾证明 $r = 1$ 时差分解仍收敛，但要求有更光滑的初值（参看注1.1及[31]）。习惯上也称 $r \leqslant 1$ 为Courant条件或CFL条件（Courant-Friedrichs-Lewy condition）。

![](images/e992f40c6326a4e992ef2351043f6cdc35daf5e9cf73b14801d9eaf1691f91c2.jpg)  
图4.3 微分方程和差分解的依存域

# 4.1.4 隐格式

为了得到恒稳定的差分格式，用第 $n - 1$ 层、 $n$ 层、 $n + 1$ 层的中心差商的权平均去逼近 $u_{xx}^{\prime \prime}$ 得到下列差分格式：

$$
\begin{array}{l} \frac {u _ {j} ^ {n + 1} - 2 u _ {j} ^ {n} + u _ {j} ^ {n - 1}}{\tau^ {2}} = a ^ {2} \left[ \theta \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} + \right. \\ \left. (1 - 2 \theta) \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} + \theta \frac {u _ {j + 1} ^ {n - 1} - 2 u _ {j} ^ {n - 1} + u _ {j - 1} ^ {n - 1}}{h ^ {2}} \right], \tag {4.23} \\ \end{array}
$$

其中 $0 \leqslant \theta \leqslant 1$ 是参数. 当 $\theta = 0$ 时就是显格式 (4.8). 实际有兴趣的参数是 $\theta = \frac{1}{4}$ , 此时差分格式可化为

$$
\left\{ \begin{array}{l} \frac {v _ {j} ^ {n + 1} - v _ {j} ^ {n}}{\tau} = a \frac {w _ {j + \frac {1}{2}} ^ {n} - w _ {j - \frac {1}{2}} ^ {n} + w _ {j + \frac {1}{2}} ^ {n + 1} - w _ {j - \frac {1}{2}} ^ {n + 1}}{2 h}, \\ \frac {w _ {j - \frac {1}{2}} ^ {n + 1} - w _ {j - \frac {1}{2}} ^ {n}}{\tau} = a \frac {v _ {j} ^ {n + 1} - v _ {j - 1} ^ {n + 1} + v _ {j} ^ {n} - v _ {j - 1} ^ {n}}{2 h}. \end{array} \right. \tag {4.24}
$$

其增长矩阵为

$$
\boldsymbol {G} (\theta) = \left[ \begin{array}{c c} \frac {1 - c ^ {2} / 4}{1 + c ^ {2} / 4} & \frac {\mathrm {i} c}{1 + c ^ {2} / 4} \\ \frac {\mathrm {i} c}{1 + c ^ {2} / 4} & \frac {1 - c ^ {2} / 4}{1 + c ^ {2} / 4} \end{array} \right], \quad c = 2 r \sin \theta .
$$

可以证明， $G(\theta)$ 的特征值按绝对值等于1,且 $\pmb{G}$ 是酉矩阵.因此 $\pmb{G}$ 的欧氏模 $\| G\| = 1$ 从而矩阵族 $\{G^n (\theta)\}$ 一致有界，故（4.22）恒稳定

# 4.1.5 数值例子

求解

$$
\begin{array}{l} u _ {t t} ^ {\prime \prime} = u _ {x x} ^ {\prime \prime}, \quad 0 <   x <   1, \quad t > 0, \\ u (0, t) = u (1, t) = 0, \quad t > 0, \\ u (x, 0) = \sin 4 \pi x, u _ {t} ^ {\prime} (x, 0) = \sin 8 \pi x, 0 <   x <   1. \\ \end{array}
$$

（精确解 $u = \cos 4\pi t\sin 4\pi x + (\sin 8\pi t\sin 8\pi x) / 8\pi)$

取空间步长 $h = \frac{1}{J}$ 时间步长 $\tau >0$ ，网比 $r = \frac{\tau}{h}$

显格式为

$$
\begin{array}{l} \frac {u _ {j} ^ {n + 1} - 2 u _ {j} ^ {n} + u _ {j} ^ {n - 1}}{\tau^ {2}} = \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}, \\ u _ {0} ^ {n} = u _ {J} ^ {n} = 0, \\ \end{array}
$$

$$
u _ {j} ^ {0} = \sin 4 \pi x _ {j}, u _ {j} ^ {1} = \sin 4 \pi x _ {j} + \tau \sin 8 \pi x _ {j}.
$$

方案I $h = \frac{1}{400} = 0.0025,\tau = \frac{1}{500} = 0.002,$ 此时 $r = \frac{4}{5}$ 计算 $t = 1,2,3,4,5$ 的差分解.

方案II $h = \tau = \frac{1}{400} = 0.0025$ ，此时 $r = 1$ 计算 $t = 1,2,3,4,5$ 的差分解

表4.1列出方案I,II的差分解的误差阶.从表中看出，方案Ⅱ的精度比方案I高很多，这是因为格式当 $r = 1$ 时截断误差有最高阶 $O\left(h^4\right),r\neq 1$ 时截断误差的阶为 $O\left(h^{2}\right)$ （参看4.1.2节的(4.7)，且时间层数受条件 $n\tau \leqslant 5$ 限制

表 4.1 方案 I 和方案 II 的差分解的误差阶  

<table><tr><td>方案</td><td>x</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td></tr><tr><td>I</td><td rowspan="2">误差阶</td><td>10-4</td><td>10-3</td><td>10-3</td><td>10-3</td><td>10-3</td></tr><tr><td>II</td><td>10-13</td><td>10-13</td><td>10-13</td><td>10-13</td><td>10-13</td></tr></table>

# 4.1.6 习题

1. 就二维波动方程导出显格式，并给出稳定性条件  
2. 证明格式 (4.24) 恒稳定  
3. 取初值 $v_{j}^{0} = (-1)^{j}, w_{j + \frac{1}{2}} = 0$ ，网比 $r = \frac{a\tau}{h} = 1$ 。求差分方程 (4.18) 的解，并用计算机验证。（答案： $v_{j}^{n} = (-1)^{n + 1}(1 - 2n), w_{j + \frac{1}{2}}^{n} = (-1)^{n + 1}2n$ ）  
4. (实习题) 利用差分格式 (4.8) — (4.10) 求下列波动方程混合边值问题的解:

$$
\left\{ \begin{array}{l l} \frac {\partial^ {2} u}{\partial t ^ {2}} - \frac {\partial^ {2} u}{\partial x ^ {2}} = 0, & 0 <   x <   1, t > 0, \\ u | _ {t = 0} = \sin \pi x, \quad \left. \frac {\partial u}{\partial t} \right| _ {t = 0} = \cos \pi x, & 0 <   x <   1, \\ u (0, t) = u (1, t) = 0, & t \geqslant 0. \end{array} \right.
$$

（精确解 $u = \sin \pi (x - t) + \sin \pi (x + t)$ ）

(1) 取 $\tau = 0.05, h = 0.1$ , 计算 $t = 0.5, 1.0, 1.5, 2.0$ 的解  
(2) 取 $\tau = h = 0.1$ , 计算 $t = 0.5, 1.0, 1.5, 2.0$ 的解

5. (实习题) 利用差分格式 (4.8) (4.9) (4.12) 重复上题的计算, 比较计算结果.

# 4.2 一阶线性双曲方程组

从本节起，我们讨论一阶线性双曲方程组的差分解法．由于构造差分格式与偏微分方程的特征及解的性质有关，所以在讨论数值解法之前，先回顾一下双曲方程组的某些基本概念（参看[4]).

# 4.2.1 双曲型方程组及其特征

设有含 $n$ 个未知函数 $\pmb{u} = (u_{1}(x,t),u_{2}(x,t),\dots ,u_{n}(x,t))$ 和 $n$ 个方程的一阶线性偏微分方程组：

$$
L _ {i} (\boldsymbol {u}) = \sum_ {j = 1} ^ {n} b _ {i j} \frac {\partial u _ {j}}{\partial t} + \sum_ {j = 1} ^ {n} a _ {i j} \frac {\partial u _ {j}}{\partial x} = c _ {i}, \quad i = 1, 2, \dots , n, \tag {4.25}
$$

其中 $b_{ij} = b_{ij}(x,t),a_{ij} = a_{ij}(x,t)$ 及 $c_{i} = c_{i}(x,t)$ 都是域 $G$ 上的光滑函数.采用矩阵、向量记号

$$
\boldsymbol {B} = \left(b _ {i j}\right) _ {n \times n}, \boldsymbol {A} = \left(a _ {i j}\right) _ {n \times n},
$$

$$
\boldsymbol {c} = \left(c _ {1}, c _ {2}, \dots , c _ {n}\right) ^ {\mathrm {T}},
$$

将(4.25)写成

$$
L (\boldsymbol {u}) = \boldsymbol {B} \frac {\partial \boldsymbol {u}}{\partial t} + \boldsymbol {A} \frac {\partial \boldsymbol {u}}{\partial x} = \boldsymbol {c}. \tag {4.26}
$$

假定矩阵 $\pmb{B}$ 有逆，则不失一般性地可设 $B = I$ （单位矩阵），只考虑如下形式的方程组：

$$
\frac {\partial \boldsymbol {u}}{\partial t} + \boldsymbol {A} \frac {\partial \boldsymbol {u}}{\partial x} = \boldsymbol {c}. \tag {4.27}
$$

定义4.1 我们说（4.27）于点 $(x,t) \in G$ 是（狭义）双曲方程组，如果矩阵 $A = A(x,t)$ 有 $n$ 个实的互异特征值

$$
\lambda_ {1} (x, t) <   \lambda_ {2} (x, t) <   \dots <   \lambda_ {n} (x, t). \tag {4.28}
$$

假若 (4.27) 于每一点 $(x, t) \in G$ 为双曲，则说它是 $G$ 上的双曲方程组。

本节总设（4.27）是双曲方程组，行向量

$$
\boldsymbol {l} ^ {(1)} (x, t), \boldsymbol {l} ^ {(2)} (x, t), \dots , \boldsymbol {l} ^ {(n)} (x, t) \tag {4.29}
$$

是矩阵 $\mathbf{A}$ 相应于 $\lambda_1, \lambda_2, \dots, \lambda_n$ 的左特征向量系，即

$$
l ^ {(i)} A = \lambda_ {i} l ^ {(i)}, \quad i = 1, 2, \dots , n.
$$

用 $l^{(i)}$ 左乘 (4.27) 两端, 得

$$
l ^ {(i)} u _ {t} ^ {\prime} + \lambda_ {i} l ^ {(i)} u _ {x} ^ {\prime} = l ^ {(i)} c, \quad i = 1, 2, \dots , n. \tag {4.30}
$$

设 $l^{(i)} = \left(l_1^{(i)}, l_2^{(i)}, \dots, l_n^{(i)}\right)$ , 则 (4.30) 即

$$
\sum_ {j = 1} ^ {n} l _ {j} ^ {(i)} \left(\frac {\partial u _ {j}}{\partial t} + \lambda_ {i} \frac {\partial u _ {j}}{\partial x}\right) = \sum_ {j = 1} ^ {n} l _ {j} ^ {(i)} c _ {j}, \quad i = 1, 2, \dots , n. \tag {4.31}
$$

在 $x, t$ 平面域 $G$ 内各点作 $n$ 个方向

$$
\mathrm {d} t: \mathrm {d} x = 1: \lambda_ {i}, \tag {4.32}
$$

或

$$
\tau_ {i}: \frac {\mathrm {d} x}{\mathrm {d} t} = \lambda_ {i}, \quad i = 1, 2, \dots , n, \tag {4.33}
$$

则沿 $\tau_{i}$

$$
\left(\frac {\mathrm {d} u _ {j}}{\mathrm {d} t}\right) _ {\tau_ {i}} = \frac {\partial u _ {j}}{\partial t} + \lambda_ {i} \frac {\partial u _ {j}}{\partial x}.
$$

于是将（4.31）化成常微分方程组：

$$
\sum_ {j = 1} ^ {n} l _ {j} ^ {(i)} \left(\frac {\mathrm {d} u _ {j}}{\mathrm {d} t}\right) _ {\tau_ {i}} = \sum_ {j = 1} ^ {n} l _ {j} ^ {(i)} c _ {j}, \quad i = 1, 2, \dots , n. \tag {4.34}
$$

由（4.33）确定的 $n$ 个不同方向称为特征方向；由特征方向确定的 $n$ 族曲线，称为特征曲线，简称特征。沿每一特征，方程组（4.27）化成常微分方程组（4.34），称为原方程的特征关系。在特征上成立的特征关系乃是利用特征概念研究双曲型方程的基础。

特征关系 (4.34) 的每一方程, 只出现未知函数沿同一方向的导数. 实际上还可进一步化简, 使每个方程只出现一个函数的方向导数. 因

$$
l _ {j} ^ {(i)} \left(\frac {\mathrm {d} u _ {j}}{\mathrm {d} t}\right) _ {\tau_ {i}} = \left(\frac {\mathrm {d} l _ {j} ^ {(i)} u _ {j}}{\mathrm {d} t}\right) _ {\tau_ {i}} - u _ {j} \left(\frac {\mathrm {d} l _ {j} ^ {(i)}}{\mathrm {d} t}\right) _ {\tau_ {i}},
$$

若令

$$
r _ {i} = r _ {i} (x, t) = \sum_ {j = 1} ^ {n} l _ {j} ^ {(i)} u _ {j}, \tag {4.35}
$$

则可将(4.34)化成

$$
\left(\frac {\mathrm {d} r _ {i}}{\mathrm {d} t}\right) _ {\tau_ {i}} = \frac {\partial r _ {i}}{\partial t} + \lambda_ {i} \frac {\partial r _ {i}}{\partial x} = \sum_ {j = 1} ^ {n} \left(l _ {j} ^ {(i)} c _ {i} + u _ {j} \left(\frac {\mathrm {d} l _ {j} ^ {(i)}}{\mathrm {d} t}\right) _ {\tau_ {i}}\right).
$$

由于（4.35）右端系数矩阵 $\pmb{L} = \left(l_j^{(i)}\right)_{n\times n}$ 有逆 $\pmb{L}^{-1} = \pmb {M} = (m_{ij})_{n\times n}$ ，故 $\pmb{u}$ 可通过 $\pmb {r} = (r_1,r_2,\dots ,r_n)^{\mathrm{T}}$ 表示为

$$
u _ {j} = \sum_ {k = 1} ^ {n} m _ {j k} r _ {k}, \quad j = 1, 2 \dots , n.
$$

于是

$$
\left(\frac {\mathrm {d} r _ {i}}{\mathrm {d} t}\right) _ {\tau_ {i}} = g _ {i} = \sum_ {j = 1} ^ {n} l _ {j} ^ {(i)} c _ {j} + \sum_ {j = 1} ^ {n} \left(\frac {\mathrm {d} l _ {j} ^ {(i)}}{\mathrm {d} t}\right) \sum_ {\tau_ {i}} ^ {n} m _ {j k} r _ {k}, i = 1, 2, \dots , n. \tag {4.36}
$$

这样就把(4.34)化成了对角方程组(4.36).新变量 $r_i(i = 1,2,\dots ,n)$ 称为Riemann不人物简介变量.

当系数矩阵 $\mathbf{A}$ 为常矩阵（因而 $l_j^{(i)}$ 与 $x,t$ 无关），且右端 $c_{i} = 0$ 时，(4.36)简化为

$$
\left(\frac {\mathrm {d} r _ {i}}{\mathrm {d} t}\right) _ {\tau_ {i}} = 0.
$$

可见 $r_i$ 沿特征

$$
\frac {\mathrm {d} x}{\mathrm {d} t} = \lambda_ {i}, \quad i = 1, 2, \dots , n
$$

是常量.

例4.1 波动方程（4.1）可化成一阶方程组：

$$
\frac {\partial \boldsymbol {u}}{\partial t} - \boldsymbol {A} \frac {\partial \boldsymbol {u}}{\partial x} = \boldsymbol {0}, \tag {4.37}
$$

$$
\boldsymbol {A} = \left[ \begin{array}{c c} 0 & a \\ a & 0 \end{array} \right],
$$

其中 $\pmb{u} = (v, w)^{\mathrm{T}}$ . $-\pmb{A}$ 的特征方程为 $\lambda^2 - a^2 = 0$ , 有二互异实根 $\pm a$ , 故 (4.37) 是双曲型方程. 特征是两族直线 $x - at = c_1, x + at = c_2$ . 特征关系为

$$
\begin{array}{l} \frac {\mathrm {d} v}{\mathrm {d} t} - \frac {\mathrm {d} w}{\mathrm {d} t} = 0, \quad \text {沿} x - a t = c _ {1}. \\ \frac {\mathrm {d} v}{\mathrm {d} t} + \frac {\mathrm {d} w}{\mathrm {d} t} = 0, \text {沿} x + a t = c _ {2}. \\ \end{array}
$$

由 (4.35) 定义的 Riemann 不变量为

$$
r _ {1} = v - w, \quad r _ {2} = v + w.
$$

于是（4.37）可写成Riemann不变量形式：

$$
\frac {\partial r _ {1}}{\partial t} + a \frac {\partial r _ {1}}{\partial x} = 0, \quad \frac {\partial r _ {2}}{\partial t} - a \frac {\partial r _ {2}}{\partial x} = 0.
$$

显然 $r_1$ 沿 $x - at = c_{1}$ 等于常数， $r_2$ 沿 $x + at = c_{2}$ 等于常数，因此

$$
r _ {1} = f (x - a t), \quad r _ {2} = g (x + a t).
$$

从而

$$
\begin{array}{l} v = \frac {1}{2} (f (x - a t) + g (x + a t)), \\ w = - \frac {1}{2} (f (x - a t) - g (x + a t)). \\ \end{array}
$$

例4.2 考虑在静止气体中小扰动(声音)传播所满足的方程组：

$$
\left\{ \begin{array}{l} \frac {\partial u}{\partial t} + \frac {c _ {0} ^ {2}}{\rho_ {0}} \frac {\partial \rho}{\partial x} = 0, \\ \frac {\partial \rho}{\partial t} + \rho_ {0} \frac {\partial u}{\partial x} = 0, \end{array} \right. \tag {4.38}
$$

其中 $u$ 和 $\rho$ 分别表示扰动后的质点速度及密度， $\rho_0$ 及 $c_{0}$ 表示静止气体的密度和音速 $(\rho_0, c_0$ 都是正常数).

现在

$$
\boldsymbol {A} = \left[ \begin{array}{c c} 0 & \frac {c _ {0} ^ {2}}{\rho_ {0}} \\ \rho_ {0} & 0 \end{array} \right],
$$

其特征值为 $\pm c_{0}$ ，所以(4.38)是双曲型方程.特征是两族直线 $x - c_0t = c_1,x + c_0t = c_2$ 特征关系为

$$
\rho_ {0} {\frac {\mathrm {d} u}{\mathrm {d} t}} + c _ {0} {\frac {\mathrm {d} \rho}{\mathrm {d} t}} = 0, \quad \text {沿} x - c _ {0} t = c _ {1},
$$

$$
\rho_ {0} \frac {\mathrm {d} u}{\mathrm {d} t} - c _ {0} \frac {\mathrm {d} \rho}{\mathrm {d} t} = 0, \quad \text {沿} x + c _ {0} t = c _ {2}.
$$

Riemann 不变量为

$$
r _ {1} = \rho_ {0} u + c _ {0} \rho , \quad r _ {2} = \rho_ {0} u - c _ {0} \rho .
$$

$r_1$ 沿 $x - c_0t = c_1$ 等于常数, $r_2$ 沿 $x + c_0t = c_2$ 等于常数

# 4.2.2 Cauchy 问题、依存域、影响域和决定域

考虑双曲方程组（4.27）的如下Cauchy问题：在线段

$$
\overline {{P _ {1} P _ {2}}} \colon t = 0, \quad a \leqslant x \leqslant b
$$

的一邻域内，求（4.27）的解 $\pmb{u}(x,t) = (u_1(x,t), u_2(x,t), \dots, u_n(x,t))^{\mathrm{T}}$ ，在 $\overline{P_1P_2}$ 上取给定的初值 $\pmb{u}^0(x)$ ，即

$$
\boldsymbol {u} (x, 0) = \boldsymbol {u} ^ {0} (x), \quad x \in \overline {{P _ {1} P _ {2}}}. \tag {4.39}
$$

由于 $x$ 轴的方向 $(\mathrm{d}x,\mathrm{d}t) = (1,0)$ 不是特征方向（因 $-\infty <  \lambda_{i} <   \infty ,i = 1,2,\dots ,n)$ ，所以Cauchy问题适定，即在 $\overline{P_1P_2}$ 一邻域内有唯一解，且解连续依赖初值(参看[3]).

现在讨论解与初值的关系. 我们称由

$$
\frac {\mathrm {d} x}{\mathrm {d} t} = \lambda_ {i} (x, t), \quad i = 1, 2, \dots , n
$$

确定的曲线族为第 $i$ 族特征.由常微分方程理论知道，同族特征不相交，所以过任一点恰有 $n$ 个不同特征经过.又 $\mathcal{X}$ 轴 $(t = 0)$ 方向不是特征方向，对任一点 $P = (x,t)(t > 0)$ 取过 $P$ 的两条特征

$$
\begin{array}{l} \tau_ {n}: \frac {\mathrm {d} x}{\mathrm {d} t} = \lambda_ {n} (x, t) = \lambda_ {\max }, \\ \tau_ {1}: \frac {\mathrm {d} x}{\mathrm {d} t} = \lambda_ {1} (x, t) = \lambda_ {\min }, \\ \end{array}
$$

假定它们与 $x$ 轴依次交于 $P_{1}$ 和 $P_{2}$ . 曲线段 $\overline{PP_1}$ 和 $\overline{PP_2}$ 与直线段 $\overline{P_1P_2}$ 围成一曲边三角形, 如图4.4(a)所示. 显然过 $P$ 的 $n$ 条特征随 $t$ 递减, 都和 $x$ 轴交于 $\overline{P_1P_2}$ . 现在从 $P$ 出发沿特征 (按 $t$ 减小方向) 积分常微分方程组 (4.36), 则知Riemann不变量 $\boldsymbol{r} = (r_1,r_2,\dots ,r_n)$ , 从而 $u = u(x,t)$ 于 $P$ 的值只和 $\overline{P_1P_2}$ 上的初值有关. 称这样的线段为点 $P$ 的依存域. 又曲边三角形 $P_{1}PP_{2}$ 围成的区域 $G$ 内任一点的 $u(x,t)$ 也由 $\overline{P_1P_2}$ 上的初值唯一决定, 故 $G$ 称为 $\overline{P_1P_2}$ 的决定域. 同时也看出, 线段 $\overline{Q_1Q_2}$ (图4.4(b))上的初值, 随 $t$ 增加可影响到上半平面过 $Q_{1}$ 的第一特征与过 $Q_{2}$ 的第 $n$ 特征之间的区域中任一点, 但不影响其他点, 故称之为 $\overline{Q_1Q_2}$ 的影响域.

![](images/0387b3ffd40daa1c32e7671163fee72e391248c34d1347b0329ca5fc632b4e1a.jpg)  
(a) 依存域

![](images/aa38688eb5ec9465edc1cc8c784cbd2ceccbeda210544ae284cc910f85cbf5f8.jpg)  
(b) 影响域  
图4.4 依存域和影响域

# 4.2.3 初边值问题

为简单计, 我们只讨论含两个未知函数 $u, v$ 的如下对角形方程组:

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} + \lambda_ {1} \frac {\partial u}{\partial x} = f _ {1}, \\ \frac {\partial v}{\partial t} + \lambda_ {2} \frac {\partial v}{\partial x} = f _ {2}, \end{array} \right. \tag {4.40}
$$

其中 $\lambda_1 < \lambda_2, f_i = a_i u + b_i v + c_i, i = 1, 2.$ 此时平面上有两族特征，我们规定特征的正向是指向 $t$ 增加的方向. 过 $(x, t)$ 的任一非特征方向 $\alpha$ 说是时向的，如果 $\alpha$ 或 $-\alpha$ 介于过此点的两正特征方向，否则就说 $\alpha$ 是空向的（参看图4.5). 在 $t > 0$ 和 $0 < x < l$ 求方程(4.40）的解 $u, v,$ 满足初值条件

$$
u (x, 0) = u ^ {0} (x), \quad v (x, 0) = v ^ {0} (x), \quad 0 \leqslant x \leqslant l \tag {4.41}
$$

和适当的边值条件, 使解存在且唯一. 如何给边值条件才是恰当的? 这里我们只作一粗略但颇有启发性的分析, 严格论证就不给了, 有兴趣的读者可参看偏微分方程的专门著作 (例如 [3]).

![](images/894aac6f026caf6632e13aa1d14c9d54bb5e34cc3a04972c6b579bf852325508.jpg)  
图4.5 时向和空向

先设方程的系数是常数，且右端 $f_{i} = 0(i = 1,2)$ ，此时两族特征都是直特征，方程(4.40）可积分，得

$$
u (x, t) = \text {常 数}, \text {当} x - \lambda_ {1} t = c _ {1} \text {时},
$$

$$
v (x, t) = \text {常 数}, \text {当} x - \lambda_ {2} t = c _ {2} \text {时}.
$$

设 $G:0 <   x <   l,t > 0.$ 分三种情形(参看图4.6):

(i) $\lambda_2 > \lambda_1 > 0$ . 此时随时间 $t$ 递增, 两族特征由左边界 $(x = 0)$ 进入 $G$ , 由右边界 $(x = l)$ 离开 $G$ (垂线 $x = 0$ 和 $x = l$ 是空向). $u, v$ 在 $G$ 的性质与左边值有关, 与右边值无关, 故 $u, v$ 的边值条件应给在左端点, 即

$$
u (0, t) = u _ {0} (t), \quad v (0, t) = v _ {0} (t).
$$

(ii) $\lambda_1 < \lambda_2 < 0$ . 与 (i) 类似 (垂线 $x = 0$ 和 $x = l$ 也是空向), $u$ 和 $v$ 的边值条件应给在右端点, 即

$$
u (l, t) = u _ {1} (t), \quad v (l, t) = v _ {1} (t).
$$

(iii) $\lambda_1 < 0 < \lambda_2$ , 此时第一族特征由右边界进入 $G$ , 第二族特征由左边界进入 $G$ (垂线 $x = 0$ 和 $x = l$ 是时向), 故边值条件应给成:

$$
u (l, t) = u _ {1} (t), \quad v (0, t) = v _ {0} (t).
$$

![](images/842d68776d98fa0e89eb502cfd6f84b5db73a1e579d610418a918e8865a8fea0.jpg)  
图4.6 三种情形

虽然上述结论是就常系数和齐右端情形给出的，但对变系数和非齐右端也成立。此外，若双曲型方程不是对角形的，则可用4.2.1小节的方法化为对角形方程组(4.36)，然后就三种不同情形对Riemann不变量 $r_i (i = 1,2)$ 给出边值条件。

# 4.2.4 习题

1. 试求下列初边值问题的解：

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} + \frac {\partial u}{\partial x} = 0, & 0 <   x <   \infty , t > 0, \\ u (x, 0) = | x - 1 |, & u (0, t) = 1. \end{array} \right.
$$

2. 试求下列初边值问题的解：

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} + (1 + x) \frac {\partial u}{\partial x} = 0, & 0 <   x <   \infty , \quad t > 0, \\ u (x, 0) = \varphi (x), & u (0, t) = 1, \quad t \geqslant 0, \end{array} \right.
$$

其中

$$
\varphi (x) = \left\{ \begin{array}{l l} 1, & 0 \leqslant x <   1, \\ 0, & x \geqslant 1. \end{array} \right.
$$

# 4.3 初值问题的差分逼近

双曲型方程与椭圆型方程、抛物型方程的一个重要区别，是双曲型方程具有特征和特征关系，其解对初值有局部依赖性质。初值的函数性质（如间断、弱间断等）将沿特征传播，因而解一般没有光滑性。在构造双曲方程的差分逼近时，应考虑这些特性。迄今已发展了许多逼近双曲型方程的差分格式，这里只介绍常见的几种，有兴趣的读者可参看文献[9, 15, 28, 31]。

# 4.3.1 迎风格式

首先考虑线性常系数方程式：

$$
\frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} = 0. \tag {4.42}
$$

这个方程虽简单，但对我们构造差分格式很有启发。我们的主要目的是构造差分格式，因此先讨论纯初值问题，然后在4.4节对初边值问题作若干注记。

沿用4.1节的记号，作(4.42)的差分逼近.按照差商代替微商的办法，自然有如下三种格式：

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j} ^ {n} - u _ {j - 1} ^ {n}}{h} = 0, \tag {4.43}
$$

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n} - u _ {j} ^ {n}}{h} = 0, \tag {4.44}
$$

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = 0. \tag {4.45}
$$

前两个方程截断误差的阶为 $O(\tau + h)$ , 第三个方程是 $O\left(\tau + h^{2}\right)$ .

从稳定性分析将会知道，这三个格式并不都可用.记

$$
r = \frac {a \tau}{h}. \tag {4.46}
$$

将(4.43)一(4.45)改写成

$$
u _ {j} ^ {n + 1} = r u _ {j - 1} ^ {n} + (1 - r) u _ {j} ^ {n}, \tag {4.47}
$$

$$
u _ {j} ^ {n + 1} = (1 + r) u _ {j} ^ {n} - r u _ {j + 1} ^ {n}, \tag {4.48}
$$

$$
u _ {j} ^ {n + 1} = u _ {j} ^ {n} + \frac {r}{2} u _ {j - 1} ^ {n} - \frac {r}{2} u _ {j + 1} ^ {n}. \tag {4.49}
$$

按Fourier方法，以 $u_{j}^{n} = v^{n}\exp (\mathrm{i}\alpha x_{j})$ $(\mathrm{i} = \sqrt{-1},x_{j} = jh,\alpha$ 是任意实参数)代到上述方程,消去公因子,分别得

$$
v ^ {n + 1} = \left(r \mathrm {e} ^ {- \mathrm {i} \alpha h} + (1 - r)\right) v ^ {n} = \lambda_ {1} v ^ {n},
$$

$$
v ^ {n + 1} = \left((1 + r) - r \mathrm {e} ^ {\mathrm {i} \alpha h}\right) v ^ {n} = \lambda_ {2} v ^ {n},
$$

$$
v ^ {n + 1} = (1 - \mathrm {i} r \sin \alpha h) v ^ {n} = \lambda_ {3} v ^ {n}.
$$

因为 $|\lambda_3| = \sqrt{1 + r^2\sin^2\alpha h}$ 对任何 $r\neq 0$ 都不满足von Neumann条件,故(4.45)恒不稳定.其次， $|\lambda_1|\leqslant 1$ 等价于 $r^2\leqslant r$ ，即

$$
\left(\frac {a \tau}{h}\right) ^ {2} \leqslant \frac {a \tau}{h},
$$

故（4.43）稳定的充要条件是

$$
a \geqslant 0, \left| \frac {a \tau}{h} \right| \leqslant 1. \tag {4.50}
$$

同理，(4.44) 稳定的充要条件是

$$
a \leqslant 0, \left| \frac {a \tau}{h} \right| \leqslant 1. \tag {4.51}
$$

这说明 $a \geqslant 0$ 时只有 (4.43) 可用, $a \leqslant 0$ 时只有 (4.44) 可用

现在我们用特征性质说明这些格式及其稳定性条件. 注意 (4.42) 的特征为

$$
\frac {\mathrm {d} x}{\mathrm {d} t} = a,
$$

特征关系是

$$
\frac {\mathrm {d} u}{\mathrm {d} t} = 0.
$$

设已知 $u_{j + 1}^n,u_j^n,u_{j - 1}^n$ ，要造出 $u_{j}^{n + 1}$ 的计算公式.如图4.7，过 $P_0(j,n + 1) = (jh,(n + 1)\tau)$ 作特征，斜率为 $\frac{\mathrm{d}t}{\mathrm{d}x} = \frac{1}{a}$ 当 $a > 0$ 时，特征偏左，与直线 $t = t_n = n\tau$ 的交点 $Q$ 位于 $Q_{0}(j,n)$ 左侧. $u$ 沿线段 $\overline{P_0Q}$ 等于常数，故 $u_{j}^{n + 1} = u_{P_{0}} = u_{Q}$ .利用 $Q_{-1},Q_{0}$ 作线性插值,得

$$
\begin{array}{l} u _ {j} ^ {n + 1} = u _ {Q} \approx \frac {u _ {Q - 1} \cdot \overline {{Q Q _ {0}}} + u _ {Q _ {0}} (h - \overline {{Q Q _ {0}}})}{h} \\ = \frac {u _ {j - 1} ^ {n} \cdot a \tau + u _ {j} ^ {n} (h - a \tau)}{h} \\ = r u _ {j - 1} ^ {n} + (1 - r) u _ {j} ^ {n}. \\ \end{array}
$$

这就是(4.47)或(4.43).稳定性条件（4.50）的第二个不等式意味着 $Q$ 应落在 $Q_{-1},Q_0$

之间, 即差分方程的依存域包含微分方程的依存域. 类似的解释也适用于 $a < 0$ , 此时过 $P_{0}$ 的特征偏右, 与 $t = t_{n}$ 的交点 $Q$ 落在 $Q_{0}$ 右侧. 利用 $Q_{0}, Q_{1}$ , 作线性插值便得到 $u_{Q}$ , 从而得出 $u_{j}^{n + 1}$ , 这就是 (4.48) 或 (4.44). 这说明差分格式 (4.43) (4.44) 与特征走向有内在联系. 用同样思想可构造变系数方程式

$$
\frac {\partial u}{\partial t} + a (x) \frac {\partial u}{\partial x} = 0
$$

的差分格式. 此时 $a$ 可能变号, 因此相应的格式为

$$
\left\{ \begin{array}{l l} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a _ {j} \frac {u _ {j} ^ {n} - u _ {j - 1} ^ {n}}{h} = 0, & a _ {j} \geqslant 0, \\ \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a _ {j} \frac {u _ {j + 1} ^ {n} - u _ {j} ^ {n}}{h} = 0, & a _ {j} <   0, \end{array} \right. \tag {4.52}
$$

其中 $a_{j} = a(x_{j})$

![](images/558b666543b9e30042589637ffb0ccbf8f1ecf6e577d2e50ed898d4d1c411e1c.jpg)

![](images/70d447f5e848f9f4cc16ed459c331f45a89ed5820a4b075e0f019dcba84dde88.jpg)  
图4.7 当 $a > 0$ 或 $a < 0$ 时的特征

这是变系数方程，为了导出稳定性条件，用局部固定系数法或视变系数为常系数法，即把 $a_{j}$ 看成与下标无关，再用Fourier方法，同样得到条件（4.50）及（4.51），只不过用 $a_{j}$ 代替那里的 $\alpha$ 就是了.

将 (4.50) (4.51) 统一写成

$$
\frac {\tau}{h} \max  _ {j} | a _ {j} | \leqslant 1. \tag {4.53}
$$

不难证明, (4.53) 是 (4.52) 稳定的充分条件. 实际上, 将 (4.52) 写成形式 (4.47) (4.48), 其中 $r = \frac{\tau}{h} a_{j}$ . 由条件 (4.53), 知 (4.47) (4.48) 右端系数非负. 当 $a_{j} \geqslant 0$ 时,

$$
\left| u _ {j} ^ {n + 1} \right| \leqslant r \left| u _ {j - 1} ^ {n} \right| + (1 - r) \left| u _ {j} ^ {n} \right| \leqslant \| U ^ {n} \| _ {\infty},
$$

当 $a_{j}\leqslant 0$ 时，

$$
\left| u _ {j} ^ {n + 1} \right| \leqslant (1 + r) \left| u _ {j} ^ {n} \right| + (- r) \left| u _ {j + 1} ^ {n} \right| \leqslant \| U ^ {n} \| _ {\infty},
$$

其中 $U^n$ 是以 $u_j^n$ 为分量的向量.总之， $\left\| U^{n + 1}\right\|_{\infty}\leqslant \left\| U^{n}\right\|_{\infty}.$ 这说明（4.52）稳定

按照气体力学的含义 $(a(x)$ 表示气流速度), 称 (4.52) 为迎风格式 (upwind scheme). 迎风格式也可推广到线性双曲方程组.

设有线性双曲型方程组：

$$
\frac {\partial \boldsymbol {u}}{\partial t} + \boldsymbol {A} \frac {\partial \boldsymbol {u}}{\partial x} = \boldsymbol {f}. \tag {4.54}
$$

设 $\lambda_1 < \lambda_2 < \dots < \lambda_m$ 是矩阵 $\mathbf{A}$ 的特征值, $l^{(1)}, l^{(2)}, \dots, l^{(m)}$ 是相应的左特征向量. 像 4.2 节那样, 用 $l^{(i)}$ 左乘 (4.54), 则得特征关系:

$$
\sum_ {k = 1} ^ {m} l _ {k} ^ {(i)} \left(\frac {\partial u _ {k}}{\partial t} + \lambda_ {i} \frac {\partial u _ {k}}{\partial x}\right) = \sum_ {k = 1} ^ {m} l _ {k} ^ {(i)} f _ {k}, \quad i = 1, 2, \dots , m. \tag {4.55}
$$

记 $\lambda_{ij} = \lambda_i(x_j)$

$$
\lambda_ {i j} \Delta^ {*} u _ {k j} ^ {n} = \left\{ \begin{array}{l l} \lambda_ {i j} \frac {u _ {k j} ^ {n} - u _ {k , j - 1} ^ {n}}{h}, & \lambda_ {i j} \geqslant 0, \\ \lambda_ {i j} \frac {u _ {k , j + 1} ^ {n} - u _ {k j} ^ {n}}{h}, & \lambda_ {i j} <   0, \end{array} \right. \tag {4.56}
$$

则逼近（4.55）的迎风差分格式是

$$
\sum_ {k = 1} ^ {m} l _ {k j} ^ {(i)} \left(\frac {u _ {k j} ^ {n + 1} - u _ {k j} ^ {n}}{\tau} + \lambda_ {i j} \Delta^ {*} u _ {k j} ^ {n}\right) = \sum_ {k = 1} ^ {m} l _ {k j} ^ {(i)} f _ {k j}, \tag {4.57}
$$

其中 $i = 1,2,\dots ,m,l_{kj}^{(i)} = l_k^{(i)}(x_j)$ 稳定性条件仍然是

$$
\left| \tau \frac {\lambda_ {i j}}{h} \right| \leqslant 1, \quad i = 1, 2, \dots , m, \quad j = 0, \pm 1, \dots . \tag {4.58}
$$

格式 (4.57) 是 Courant, Isaacson 和 Rees (1952) 提出的，并证明了在条件 (4.58) 下差分解对光滑解的收敛性。

# 4.3.2 积分守恒差分格式

前面介绍的迎风格式是根据特征走向构造的向前或向后差分格式. 现在从积分守恒方程出发构造差分格式

所谓守恒方程是指如下散度型偏微分方程：

$$
\frac {\partial u}{\partial t} + \frac {\partial f (x , u)}{\partial x} = 0. \tag {4.59}
$$

设 $G$ 是 $xt$ 平面任一有界域, 据 Green 公式,

$$
\iint_ {G} \left(\frac {\partial u}{\partial t} + \frac {\partial f}{\partial x}\right) \mathrm {d} x \mathrm {d} t = \int_ {\Gamma} (f \mathrm {d} t - u \mathrm {d} x),
$$

其中 $\Gamma = \partial G$ (取逆时针方向). 于是可将 (4.59) 写成积分守恒形式

$$
\int_ {\Gamma} (f \mathrm {d} t - u \mathrm {d} x) = 0. \tag {4.60}
$$

我们先从（4.60）出发构造熟知的Lax-Friedrichs格式.设网格如图4.8，取 $G$ 为以 $A(j + 1,n),B(j + 1,n + 1),C(j - 1,n + 1)$ 和 $D(j - 1,n)$ 为顶点的开矩形， $\Gamma = \overline{ABCDA}$ 为其边界(取逆时针方向)，则

$$
\int_ {\Gamma} (f \mathrm {d} t - u \mathrm {d} x) = \int_ {D A} (- u) \mathrm {d} x + \int_ {B C} (- u) \mathrm {d} x + \int_ {A B} f \mathrm {d} t + \int_ {C D} f \mathrm {d} t. \tag {4.61}
$$

右端第一个积分用梯形公式，第二个积分用中矩形公式，第三、四两个积分用下矩形公式，则由（4.60）（4.61）得Lax-Friedrichs格式：

$$
\frac {u _ {j} ^ {n + 1} - \frac {1}{2} \left(u _ {j - 1} ^ {n} + u _ {j + 1} ^ {n}\right)}{\tau} + \frac {f _ {j + 1} ^ {n} - f _ {j - 1} ^ {n}}{2 h} = 0, \tag {4.62}
$$

其中 $f_{j}^{n} = f\left(x_{j},u_{j}^{n}\right)$ .Lax-Friedrichs格式截断误差的阶是 $O\left(\tau +h^{2}\right)$ .特别地，当 $f =$ au时,Lax-Friedrichs格式相当于(4.45）中的 $u_{j}^{n}$ 代以 $\frac{1}{2}\left(u_{j - 1}^n +u_{j + 1}^n\right)$ 的结果.我们知道,(4.45)恒不稳定,但由Fourier方法可知,Lax-Friedrichs格式稳定的充要条件是

$$
\frac {| a | \tau}{h} \leqslant 1. \tag {4.63}
$$

![](images/094b467a73a7f639ab0a0ce69059209d803e9a2e1f967f239807b34a0e0fc9ce.jpg)  
图4.8 网格构造

若在上述推导过程中将 $u$ 换成向量 $\mathbf{U} = (u_{1}, u_{2}, \dots, u_{m})^{\mathrm{T}}$ , $f$ 换成向量 $\mathbf{f} = (f_{1}, f_{2}, \dots, f_{m})^{\mathrm{T}}$ , 则得到方程组的 Lax-Friedrichs 格式.

现在由积分守恒方程导出另一种所谓盒式格式 (box scheme). 如图 4.9, 取 $G$ 为以网点 $A(j,n), B(j,n+1), C(j-1,n+1)$ 和 $D(j-1,n)$ 为顶点的矩形, $\Gamma$ 是 $G$ 的边界. 此时积分型方程仍具有形式 (4.61). 右端各项积分用梯形公式近似, 则得

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + \frac {u _ {j - 1} ^ {n + 1} - u _ {j - 1} ^ {n}}{\tau} + \frac {f _ {j} ^ {n} - f _ {j - 1} ^ {n}}{h} + \frac {f _ {j} ^ {n + 1} - f _ {j - 1} ^ {n + 1}}{h} = 0. \tag {4.64}
$$

特别地, 以 $f = au$ 代入, 并令 $r = \frac{\tau}{h}$ , 得

$$
(1 + a u) u _ {j} ^ {n + 1} + (1 - a r) u _ {j - 1} ^ {n + 1} = (1 - a r) u _ {j} ^ {n} + (1 + a r) u _ {j - 1} ^ {n}. \tag {4.65}
$$

当 $a > 0$ 时，边值给在左端，计算由左向右进行；当 $a < 0$ 时，边值给在右端，计算由右向左进行，这些计算都是显的。用Fourier方法还可证明（4.65）恒稳定。若 $a$ 依赖 $x$ ，则

用 $a_{j - \frac{1}{2}}$ 代替（4.65）中的 $a$ . 像迎风格式一样，盒式格式也可用于双曲方程组．此外，若视（4.64）是在 $\left(x_{j - \frac{1}{2}}, t_{n + \frac{1}{2}}\right)$ 的差分逼近，则截断误差的阶为 $O\left(\tau^{2} + h^{2}\right)$

![](images/dfb41424061e49acb139cc5b89143f8e7cbf267333779be0945284bf7dda1997.jpg)  
图4.9 盒式格式

# 4.3.3 粘性差分格式

粘性差分格式的构造分两步. 先在双曲方程中引进一带二阶空间导数的小参数项, 称为粘性项, 使之成为一带小参数的抛物型方程, 例如

$$
\frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} = \varepsilon \frac {\partial^ {2} u}{\partial x ^ {2}} \quad (\varepsilon > 0); \tag {4.66}
$$

然后构造逼近相应抛物方程的差分格式. 自然要求 $\varepsilon \rightarrow 0$ , 当 $\tau \rightarrow 0$ 时

许多逼近双曲方程的差分格式可看作粘性差分格式.例如迎风格式（4.52）可改写为

$$
\begin{array}{l} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a _ {j} \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = \frac {h}{2} a _ {j} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}, \quad a _ {j} \geqslant 0, \\ \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a _ {j} \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = - \frac {h}{2} a _ {j} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}, \quad a _ {j} <   0. \\ \end{array}
$$

或写成统一形式

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a _ {j} \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = \frac {h}{2} | a _ {j} | \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}. \tag {4.67}
$$

这相当于下列带小参数的抛物方程的中心差分逼近：

$$
\frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} = \frac {h}{2} | a | \frac {\partial^ {2} u}{\partial x ^ {2}}. \tag {4.68}
$$

对于线性双曲方程组 (4.54), 设 $A = S\Lambda S^{-1}, \Lambda = \mathrm{diag}(\lambda_1, \lambda_2, \dots, \lambda_n)$ 是由 $A$ 的特征值组成的对角矩阵, 则可将逼近 (4.54) 的迎风格式写成:

$$
\frac {\boldsymbol {u} _ {j} ^ {n + 1} - \boldsymbol {u} _ {j} ^ {n}}{\tau} + \boldsymbol {A} \frac {\boldsymbol {u} _ {j + 1} ^ {n} - \boldsymbol {u} _ {j - 1} ^ {n}}{2 h} = \frac {h}{2} \boldsymbol {S} | \boldsymbol {\Lambda} | \boldsymbol {S} ^ {- 1} \frac {\boldsymbol {u} _ {j + 1} ^ {n} - 2 \boldsymbol {u} _ {j} ^ {n} + \boldsymbol {u} _ {j - 1} ^ {n}}{h ^ {2}} \tag {4.69}
$$

其中 $|\pmb {\Lambda}| = \mathrm{diag}(|\lambda_1|,|\lambda_2|,\dots ,|\lambda_n|).$

Lax-Friedrichs格式（4.62）也可改写成：

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + \frac {f _ {j + 1} ^ {n} - f _ {j - 1} ^ {n}}{2 h} = \frac {h ^ {2}}{2 \tau} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}. \tag {4.70}
$$

可看成是带小参数的抛物方程：

$$
\frac {\partial u}{\partial t} + \frac {\partial f}{\partial x} = \frac {h}{2 r} \frac {\partial^ {2} u}{\partial x ^ {2}} \tag {4.71}
$$

中心差分化的结果，网比 $r = \frac{\tau}{h}$ 固定

比较 (4.68) 和 (4.71) 可知, 小参数的不同取法可导出各种不同的差分格式, 所以选取小参数是构造粘性差分格式的关键. 对于实际问题, 可以直接取自然粘性项, 也可以根据物理上的某些考虑构造人工粘性项. von Neumann 和 Richtmyer 引进的粘性项就属于后一种 (参看 [28] §12.10). 现在以 Lax-Wendroff 格式为例, 介绍另一种引进人工粘性项的方法.

为简便起见, 不妨设 (4.59) 中的 $f = f(u)$ (不显含 $x$ ). 将 $u$ 关于时间变量 $t$ 展开, 有

$$
u \left(x _ {j}, t _ {n + 1}\right) = u \left(x _ {j}, t _ {n}\right) + \tau \left(\frac {\partial u}{\partial t}\right) _ {j} ^ {n} + \frac {\tau^ {2}}{2} \left(\frac {\partial^ {2} u}{\partial t ^ {2}}\right) _ {j} ^ {n} + \dots .
$$

利用 (4.59) 将关于 $t$ 的偏导数换成关于 $x$ 的偏导数，并注意

$$
\frac {\partial^ {2} u}{\partial t ^ {2}} = - \frac {\partial}{\partial t} \frac {\partial f}{\partial x} = - \frac {\partial}{\partial x} \frac {\partial f}{\partial t} = - \frac {\partial}{\partial x} \left(f ^ {\prime} (u) \frac {\partial u}{\partial t}\right) = \frac {\partial}{\partial x} \left(f ^ {\prime} (u) \frac {\partial f}{\partial x}\right),
$$

则

$$
u \left(x _ {j}, t _ {n + 1}\right) = u \left(x _ {j}, t _ {n}\right) - \tau \left(\frac {\partial f}{\partial x}\right) _ {j} ^ {n} + \frac {\tau^ {2}}{2} \left(\frac {\partial}{\partial x}\right) \left(f ^ {\prime} (u) \frac {\partial f}{\partial x}\right) _ {j} ^ {n} + \dots .
$$

然后略去余项，并用中心差商代替对 $x$ 的偏导数，则得Lax-Wendroff格式：

$$
u _ {j} ^ {n + 1} = u _ {j} ^ {n} - \frac {1}{2} \frac {\tau}{h} \left(f _ {j + 1} ^ {n} - f _ {j - 1} ^ {n}\right) + \frac {1}{2} \left(\frac {\tau}{h}\right) ^ {2} \left[ a _ {j + \frac {1}{2}} ^ {n} \left(f _ {j + 1} ^ {n} - f _ {j} ^ {n}\right) - a _ {j - \frac {1}{2}} ^ {n} \left(f _ {j} ^ {n} - f _ {j - 1} ^ {n}\right) \right], \tag {4.72}
$$

其中 $a_{j + \frac{1}{2}}^n = f'\left(\frac{1}{2} u_{j + 1}^n +\frac{1}{2} u_j^n\right)$ 特别地，当 $a$ 是常数时，(4.72）为

$$
u _ {j} ^ {n + 1} = u _ {j} ^ {n} - \frac {1}{2} a \frac {\tau}{h} \left(u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}\right) + \frac {1}{2} \left(a \frac {\tau}{h}\right) ^ {2} \left(u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}\right). \tag {4.73}
$$

Lax-Wendroff格式可看成是带粘性项方程

$$
\frac {\partial u}{\partial t} + \frac {\partial f}{\partial x} = \frac {\tau}{2} \frac {\partial}{\partial x} \left(f ^ {\prime} (u) \frac {\partial f}{\partial x}\right) \tag {4.74}
$$

的中心差分格式. 显然, 格式 (4.72) 的截断误差的阶是 $O\left(\tau^{2} + h^{2}\right)$ . 由Fourier方法可知 (4.73) 的增长因子是

$$
G = 1 - \operatorname {i r a} \sin (\alpha h) - (r a) ^ {2} (1 - \cos (\alpha h)),
$$

其中 $r = \frac{\tau}{h}$ 为使 $|G|\leqslant 1$ ，必须且只需 $r|a|\leqslant 1.$ 这就是稳定性条件

为了避免出现导数 $a = f'(u)$ ，在网格中心引进过渡值

$$
u _ {j + \frac {1}{2}} ^ {n + \frac {1}{2}} = \frac {1}{2} \left(u _ {j + 1} ^ {n} + u _ {j} ^ {n}\right) - \frac {r}{2} \left(f _ {j + 1} ^ {n} - f _ {j} ^ {n}\right), \tag {4.75}
$$

再由下式得到最终值：

$$
u _ {j} ^ {n + 1} = u _ {j} ^ {n} - r \left(f _ {j + \frac {1}{2}} ^ {n + \frac {1}{2}} - f _ {j - \frac {1}{2}} ^ {n + \frac {1}{2}}\right). \tag {4.76}
$$

称（4.75）和（4.76）为两步Lax-Wendroff法，它仍然有二阶截断误差，且当 $a = f'(x) =$ 常数时可转化成(4.73).

粘性差分格式(4.70)和(4.72)也可直接推广到双曲方程组，只需将那里的 $u$ 和 $f$ 换成 $m$ 维向量函数

# 4.3.4 其他差分格式

下面列出逼近（4.42）的其他一些差分格式及相关结果，这些格式各有其特点。以下用 $r = \frac{a\tau}{h}$ 表示网比，并令

$$
\Delta_ {+} u _ {j} = u _ {j + 1} - u _ {j}, \quad \Delta_ {-} u _ {j} = u _ {j} - u _ {j - 1}, \quad \Delta_ {0} u _ {j} = u _ {j + 1} - u _ {j - 1}.
$$

1. Beam-Warming 格式:

$$
\left\{ \begin{array}{l} u _ {j} ^ {*} = u _ {j} ^ {n} - r \Delta_ {-} u _ {j} ^ {n}, \\ u _ {j} ^ {n + 1} = \frac {1}{2} \left(u _ {j} ^ {n} + u _ {j} ^ {*} - r \Delta_ {-} u _ {j} ^ {*} - r \Delta_ {-} \Delta_ {+} u _ {j - 1} ^ {n}\right). \end{array} \right. \tag {4.77}
$$

稳定性条件： $0\leqslant |r|\leqslant 2$ ，截断误差阶： $O(\tau^2) + O(\tau h) + O(h^2)$

2. MacCormack 格式:

$$
\left\{ \begin{array}{l} u _ {j} ^ {*} = u _ {j} ^ {n} - r \Delta_ {+} u _ {j} ^ {n}, \\ u _ {j} ^ {n + 1} = \frac {1}{2} \left(u _ {j} ^ {n} + u _ {j} ^ {*} - r \Delta_ {-} u _ {j} ^ {*}\right). \end{array} \right. \tag {4.78}
$$

稳定性条件： $|r|\leqslant 1$ ，截断误差阶： $O\left(\tau^{2}\right) + O\left(h^{2}\right)$

3. 隐式迎风格式：

$$
\left\{ \begin{array}{l} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j} ^ {n + 1} - u _ {j - 1} ^ {n + 1}}{h} = 0, \quad a \geqslant 0, \\ \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n + 1} - u _ {j} ^ {n + 1}}{h} = 0, \quad a <   0 \end{array} \right. \tag {4.79}
$$

恒稳定. 截断误差阶: $O(\tau) + O(h)$ .

4. 隐式中心格式：

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n + 1} - u _ {j - 1} ^ {n + 1}}{2 h} = 0 \tag {4.80}
$$

恒稳定. 截断误差阶: $O(\tau) + O(h^{2})$ .

5. 跳蛙 (leap-frog) 格式:

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n - 1}}{2 \tau} + a \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = 0. \tag {4.81}
$$

稳定条件: $\left|\frac{a\tau}{h}\right| \leqslant 1$ , 截断误差阶: $O\left(\tau^{2}\right) + O\left(h^{2}\right)$ .

# 4.3.5 习题

1. 逼近方程

$$
\frac {\partial u}{\partial t} + \frac {\partial (a u)}{\partial x} = 0 \quad (a = a (x))
$$

的另一形式的迎风格式为

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + \frac {a _ {j + \frac {1}{2}} u _ {j + \frac {1}{2}} ^ {n} - a _ {j - \frac {1}{2}} u _ {j - \frac {1}{2}} ^ {n}}{h} = 0, \tag {4.82}
$$

其中

$$
u _ {j + \frac {1}{2}} ^ {n} = \left\{ \begin{array}{l l} u _ {j + 1} ^ {n}, & a _ {j + \frac {1}{2}} <   0, \\ u _ {j} ^ {n}, & a _ {j + \frac {1}{2}} > 0, \end{array} \right. \tag {4.83}
$$

$$
u _ {j - \frac {1}{2}} ^ {n} = \left\{ \begin{array}{l l} u _ {j} ^ {n}, & a _ {j - \frac {1}{2}} <   0, \\ u _ {j - 1} ^ {n}, & a _ {j - \frac {1}{2}} > 0. \end{array} \right. \tag {4.84}
$$

因 $a$ 表示流速，所以也称(4.82)—(4.84)为偏上游格式(upstream scheme).证明当 $\tau \leqslant \frac{h}{\sup |a(x)|}$ 时格式稳定

2. 利用特征构造逼近方程 (4.42) 的两种隐式迎风格式. 第一种, 如图 4.10(a), 过 $P_{0}$ 作特征, 若特征与对角线 $\overline{Q_{0}P_{-1}}$ 相交, 则利用 $u_{Q_{0}}, u_{P_{-1}}$ 作交点的线性插值, 并取 $u_{P_{0}}$ 等于交点值. 若特征与对角线 $\overline{Q_{0}P_{1}}$ 相交, 可类似地作 $u_{P_{0}}$ . 第二种, 如图 4.10(b), 过 $P_{0}$ 作特征, 则特征必与线段 $\overline{Q_{0}Q_{1}}, \overline{Q_{1}P_{1}}, \overline{Q_{0}Q_{-1}}$ 和 $\overline{Q_{-1}P_{-1}}$ 之一相交, 例如和 $\overline{Q_{1}P_{1}}$ 相交, 则利用 $u_{Q_{1}}, u_{P_{1}}$ 作交点的线性插值, 并取 $u_{P_{0}}$ 等于交点值, 余类推. 试导出计算公式并研究稳定性.

![](images/9d26e403fb633f9b52751ad96cec3b7663621a1d6d4db9f6f8c2ce76864337de.jpg)  
(a)第一种

![](images/6c7834893f097b15e2d4e7ce8899a88377e071332f0ca2e11282c82c377f4181.jpg)  
(b) 第二种  
图4.10 两种隐式迎风格式

3. (实习题) 利用显和隐的迎风格式计算：

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} - 2 \frac {\partial u}{\partial x} = 0, & x \in (0, 1), \quad t > 0, \\ u (x, 0) = 1 + \sin 2 \pi x, & x \in [ 0, 1 ], \\ u (1, t) = 1. \end{array} \right.
$$

取 $h = 0.1, \tau = 0.02$ ，算出解在 $t = 0.1, 0.5$ 的值

4. (实习题) 求解初边值问题:

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} + \frac {\partial u}{\partial x} = 0, & 0 <   x <   \infty , \quad t > 0, \\ u (x, 0) = | x - 1 |, & u (0, t) = 1. \end{array} \right. \tag {4.85}
$$

(1) 利用格式 (4.43), 取 $\tau = h = 0.5$ , 计算 $t = 1,2,3,4,5$ 的解.  
(2) 利用格式 (4.79), 取 $\tau = 1, h = 0.5$ , 计算 $t = 1, 2, 3, 4, 5$ 的解.

注 方程 (4.85) 的特征族为 $x - t = c$ . 据初边值条件, 精确解为

$$
u (x, t) = \left\{ \begin{array}{l l} 1, & - \infty <   x - t <   0, \\ 1 - (x - t), & 0 <   x - t <   1, \\ (x - t) - 1, & 1 <   x - t <   \infty . \end{array} \right.
$$

# *4.4 初边值问题和对流占优扩散方程

# 4.4.1 初边值问题

以模型问题

$$
\frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} = 0, \quad x \in (0, 1), \quad t > 0, \tag {4.86}
$$

$$
u (x, 0) = \varphi (x), \quad x \in [ 0, 1 ] \tag {4.87}
$$

为例介绍边值条件几种给法及其逼近方法

1. 周期边值条件 是指边值条件：

$$
u (0, t) = u (1, t).
$$

设初值函数也以1为周期： $\varphi (0) = \varphi (1)$ .则可将此初边值问题以周期1扩展到 $x$ 轴，使其成为纯初值问题.设 $a > 0,r = a\frac{\tau}{h}$ ，则求解它的迎风格式为

$$
\begin{array}{l} u _ {j} ^ {n + 1} = u _ {j} ^ {n} - r \left(u _ {j} ^ {n} - u _ {j - 1} ^ {n}\right), \quad j = 1, 2, \dots , J, (4.88) \\ u _ {0} ^ {n + 1} = u _ {J} ^ {n + 1}, \\ u _ {j} ^ {0} = \varphi_ {j} = \varphi (x _ {j}), \quad j = 1, 2, \dots , J, (4.89) \\ n = 0, 1, \dots . \\ \end{array}
$$

稳定性条件为 $r \leqslant 1$

若用Lax-Friedrichs格式，则

$$
\begin{array}{l} u _ {j} ^ {n + 1} = \frac {1}{2} (1 + r) u _ {j + 1} ^ {n} + \frac {1}{2} (1 - r) u _ {j - 1} ^ {n}, \quad j = 1, 2, \dots , J - 1, \\ u _ {0} ^ {n + 1} = u _ {J} ^ {n + 1}, \quad n = 0, 1, \dots , \tag {4.90} \\ u _ {j} ^ {0} = \varphi_ {j}, \quad j = 1, 2, \dots , J - 1. \\ \end{array}
$$

因 $u_{-1}^{n} = u_{J - 1}^{n}$ ，所以

$$
u _ {0} ^ {n + 1} = \frac {1}{2} (1 + r) u _ {1} ^ {n} + \frac {1}{2} (1 - r) u _ {J - 1} ^ {n}.
$$

2. Dirichlet 条件 此时要按照 $a$ 的符号配置边值。根据特征的走向，边值应如下配置：

$$
\begin{array}{l} u (0, t) = u _ {0} (x), \quad a > 0 \\ u (1, t) = u _ {1} (x), \quad a <   0. \\ \end{array}
$$

迎风格式应为

$$
\begin{array}{l} u _ {j} ^ {n + 1} = u _ {j} ^ {n} - r \left(u _ {j} ^ {n} - u _ {j - 1} ^ {n}\right), \quad a > 0, \\ u _ {j} ^ {n + 1} = u _ {j} ^ {n} - r \left(u _ {j + 1} ^ {n} - u _ {j} ^ {n}\right), \quad a <   0. \\ \end{array}
$$

另一个合理的选择是采用隐式迎风格式（见(4.79)）：

$$
\left\{ \begin{array}{l} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j} ^ {n + 1} - u _ {j - 1} ^ {n + 1}}{h} = 0, \quad a > 0, \\ \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n + 1} - u _ {j} ^ {n + 1}}{h} = 0, \quad a <   0. \end{array} \right.
$$

此格式恒稳定，截断误差的阶为 $O(\tau + h)$ 。虽是隐格式，但可显式求解。

3. 数值边界条件 若用空间中心差分格式逼近 (4.86) 的 Dirichlet 边值问题, 例如用 Lax-Friedrichs 格式或 Lax-Wendrof 格式, 就会发现还缺少一个边值条件, 因此需适当补充一个条件, 称为数值边界条件. 例如用 Lax-Wendrof 格式

$$
u _ {j} ^ {n + 1} = u _ {j} ^ {n} - \frac {1}{2} r \left(u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}\right) + \frac {1}{2} r ^ {2} \left(u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}\right), \tag {4.91}
$$

逼近 (4.86), 其中 $a < 0, x \in (0,1)$ , 边值条件为 $u(1,t) = 1$ . 此时应在 $x = 0$ 补充给一边值条件. 通常可按 (偏右) 迎风格式给为

$$
u _ {0} ^ {n + 1} = u _ {0} ^ {n} - r \left(u _ {1} ^ {n} - u _ {0} ^ {n}\right), \tag {4.92}
$$

或

$$
u _ {0} ^ {n + 1} = u _ {0} ^ {n} - r \left(u _ {1} ^ {n + 1} - u _ {0} ^ {n + 1}\right).
$$

也可用外推数值边值条件：

$$
u _ {0} ^ {n + 1} - u _ {1} ^ {n + 1} = 0, \tag {4.93}
$$

或

$$
u _ {0} ^ {n + 1} - 2 u _ {1} ^ {n + 1} + u _ {2} ^ {n + 1} = 0 \tag {4.94}
$$

等, 但不能任意指定 $u_0^n$ 值 (参看 [31] 的第 8 章).

# 4.4.2 对流占优扩散方程

所谓对流占优扩散方程，是指带对流项的抛物方程：

$$
\frac {\partial u}{\partial t} + b \frac {\partial u}{\partial x} = a \frac {\partial^ {2} u}{\partial x ^ {2}}, \tag {4.95}
$$

其中 $a, b$ 是常数， $a > 0$ ，而 $a \ll |b|$ （即 $|b|$ 相对 $a$ 充分大）。此时 (4.95) 虽是抛物方程，但其解却具双曲性质，构造差分格式时需考虑解的这一性质。先处理方程 (4.95) 的左端。令

$$
\alpha = \left(1 + b ^ {2}\right) ^ {\frac {1}{2}}.
$$

与 $u_{t} + bu_{x}$ 相伴的特征方向为

$$
\boldsymbol {\nu} = \left(\frac {1}{\alpha}, \frac {b}{\alpha}\right),
$$

沿 $\pmb{\nu}$ 的方向导数

$$
\frac {\partial}{\partial \nu} = \frac {1}{\alpha} \frac {\partial}{\partial t} + \frac {b}{\alpha} \frac {\partial}{\partial x}.
$$

所以（4.95）可写成形式：

$$
\alpha \frac {\partial u}{\partial \nu} = a \frac {\partial^ {2} u}{\partial x ^ {2}}. \tag {4.96}
$$

取时间步长 $\tau > 0$ ，沿 $t$ 轴取节点 $t = t_n = n\tau, n = 1,2,\dots$ 。由 $(x,t_{n})$ 出发的特征（方向为 $\pmb{\nu}$ ）与直线 $t = t_{n - 1}$ 交于

$$
\bar {x} = x - b \tau ,
$$

参看图4.11. 自然用下式逼近沿特征方向的导数：

$$
\alpha \frac {\partial u}{\partial \nu} \approx \alpha \frac {u (x , t _ {n}) - u (\bar {x} , t _ {n - 1})}{[ (x - \bar {x}) ^ {2} + \tau^ {2} ] ^ {\frac {1}{2}}} = \frac {u (x , t _ {n}) - u (\bar {x} , t _ {n - 1})}{\tau}.
$$

于是可用如下方程近似替代 (4.96)

$$
\frac {u \left(x , t _ {n}\right) - u \left(\bar {x} , t _ {n - 1}\right)}{\tau} = a \frac {\partial^ {2} u}{\partial x ^ {2}}. \tag {4.97}
$$

![](images/892eb46099147c858dff527d7b533cec492a4290ca9f657f8f6e787f54079fa5.jpg)  
图4.11 特征线

取 $x = x_{j},\tau >0$ 充分小，使 $(\bar{x},t_{n - 1})$ 位于 $(x_{j - 1},t_{n - 1}),(x_{j + 1},t_{n - 1})$ 之间，则 $b\geqslant 0$ 时位于 $(x_{j},t_{n - 1})$ 左侧， $b <   0$ 时位于 $(x_{j},t_{n - 1})$ 右侧（见图4.11).今设 $b\geqslant 0,$ 与4.3.1小节类似，在 $(x_{j - 1},t_{n - 1}),(x_j,t_{n - 1})$ 之间以 $u_{j - 1}^{n - 1},u_j^{n - 1}$ 为型值作线性插值，得

$$
\begin{array}{l} u \left(\bar {x}, t _ {n - 1}\right) = \frac {x _ {j} - \bar {x}}{h} u _ {j - 1} ^ {n - 1} + \frac {\bar {x} - x _ {j - 1}}{h} u _ {j} ^ {n - 1} \tag {4.98} \\ = \frac {b \tau}{h} u _ {j - 1} ^ {n - 1} + \frac {h - b \tau}{h} u _ {j} ^ {n - 1}. \\ \end{array}
$$

以之代到 (4.97) 左端, 并用二阶中心差商代替右端的二阶偏微商, 就得到逼近 (4.95) 的迎风差分格式

$$
\frac {u _ {j} ^ {n} - u _ {j} ^ {n - 1}}{\tau} + b \frac {u _ {j} ^ {n - 1} - u _ {j - 1} ^ {n - 1}}{h} = a \frac {u _ {j + 1} ^ {n - 1} - 2 u _ {j} ^ {n - 1} + u _ {j - 1} ^ {n - 1}}{h ^ {2}}. \tag {4.99}
$$

同样，当 $b < 0$ 时有逼近（4.95）的迎风差分格式：

$$
\frac {u _ {j} ^ {n} - u _ {j} ^ {n - 1}}{\tau} + b \frac {u _ {j + 1} ^ {n - 1} - u _ {j} ^ {n - 1}}{h} = a \frac {u _ {j + 1} ^ {n - 1} - 2 u _ {j} ^ {n - 1} + u _ {j - 1} ^ {n - 1}}{h ^ {2}}. \tag {4.100}
$$

若 $a > |b|$ ，则用中心差分格式：

$$
\frac {u _ {j} ^ {n} - u _ {j} ^ {n - 1}}{\tau} + b \frac {u _ {j + 1} ^ {n - 1} - u _ {j - 1} ^ {n - 1}}{h} = a \frac {u _ {j + 1} ^ {n - 1} - 2 u _ {j} ^ {n - 1} + u _ {j - 1} ^ {n - 1}}{h ^ {2}}. \tag {4.101}
$$

若 $|b|\gg a$ 时仍用(4.101)，则可能出现不应有的振荡.利用Fourier方法等技巧可以证明，(4.101）稳定的充要条件是（见[1])

$$
\frac {1}{2} \left(\frac {b \tau}{h}\right) ^ {2} \leqslant \frac {a \tau}{h ^ {2}} \leqslant \frac {1}{2}, \tag {4.102}
$$

或

$$
\frac {1}{2} r _ {1} ^ {2} \leqslant r \leqslant \frac {1}{2},
$$

其中

$$
r _ {1} = \frac {b \tau}{h}, r = \frac {a \tau}{h ^ {2}}.
$$

由(4.102)得 $\tau \leqslant \frac{2a}{|b|}$ 若 $|b|\gg a$ ，这条件很难满足

差分格式 (4.99) (4.100) 的截断误差的阶为 $O(\tau + h)$ . 如果用 $(x_{j-1}, t_{n-1}), (x_j, t_{n-1})$ , $(x_{j+1}, t_{n-1})$ 和型值 $u_{j-1}^{n-1}, u_j^{n-1}, u_{j+1}^{n-1}$ 的二次插值代替线性插值 (4.98)，就可得到截断误差为 $O(\tau + h^2)$ 的格式.

# 4.4.3 数值例子

求解对流占优扩散方程（利用Fourier方法等技巧）

$$
u _ {t} ^ {\prime} + b u _ {x} ^ {\prime} = a u _ {x x} ^ {\prime \prime}, \quad t > 0, \quad 0 <   x <   \infty , \tag {4.103}
$$

其中 $a, b > 0$ 是常数. 定解条件为

$$
u (x, 0) = 0, \quad x > 0, \tag {4.104}
$$

$$
u (0, t) = u _ {0}, \quad t \geqslant 0, \quad u _ {0} \text {为 正 常 数}, \tag {4.105}
$$

$$
u (\infty , t) = 0, \quad t \geqslant 0. \tag {4.106}
$$

精确解为

$$
u (x, t) = \frac {u _ {0}}{2} \left\{\operatorname {e r f c} \left(\frac {x - b t}{2 \sqrt {a t}}\right) + \exp \left(\frac {b x}{a}\right) \operatorname {e r f c} \left(\frac {x + b t}{2 \sqrt {a t}}\right) \right\}, \operatorname {e r f c} (x) = \frac {2}{\sqrt {\pi}} \int_ {x} ^ {\infty} \mathrm {e} ^ {- t ^ {2}} \mathrm {d} t.
$$

现在采用隐式迎风格式 (参看 (4.99)):

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + b \frac {u _ {j} ^ {n + 1} - u _ {j - 1} ^ {n + 1}}{h} = a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}}.
$$

(i) 设 $a = 1, b = 1$ . 取 $h = 0.5, \tau = 0.25$ . 此时 $r = \frac{\tau}{h^2} = 1$ . 记 $t_n = n\tau$ , 计算 $t = t_1, t_2, \dots, t_{50}$ 的差分解, 没有发生振荡, $t = t_{50}$ 的曲线如图 4.12, 实线为精确解, 虚线为近似解.

![](images/831a0034e85ac3cb78874f7631ba16537df89c68b230f31bf04125fd61126116.jpg)  
图4.12 差分解曲线

(ii) 设 $a = 1, b = 10$ . 取 $h = 0.5, \tau = 0.25$ . 此时 $r = 1$ . 记 $t_n = n\tau$ . 计算 $t = t_1, t_2, \dots, t_{50}$ 的差分解, 没有发生振荡, $t = t_{50}$ 的曲线如图 4.13, 实线为精确解, 虚线为近似解. 可见迎风格式计算含大对流项方程仍稳定. 我们曾用稳定的显中心差分格式求解问题 (4.103) — (4.106), 当对流项 $b = 10$ 时解出现不应有的振荡.

![](images/1c312a84eccd6f69205d56b3527673a729d469d3b3afb47d6055e138e059df4f.jpg)  
图4.13 差分解曲线

# 4.4.4 习题

1. (实习题) 用迎风格式计算求解初边值问题

$$
\begin{array}{l} u _ {t} ^ {\prime} - u _ {x} ^ {\prime} = 0, \quad x \in (0, 1), \\ u (x, 0) = \sin^ {4 0} \pi x, \quad x \in [ 0, 1 ], \\ u (0, t) = u (1, t), \quad t \geqslant 0. \\ \end{array}
$$

此问题的精确解为 $u(x,t) = \sin^{40}\pi (x + t)$

(a) 取 $h = 0.05, \tau = 0.04$ , 算出 $t = 0.00, 0.12, 0.20, 0.80$ 的值.  
(b) 画出 (a) 的图形, 观察峰值位置的变化, 与精确解峰值位置比较.

2. (实习题) 用Lax-Wendroff格式求解方程

$$
\begin{array}{l} u _ {t} ^ {\prime} - 2 u _ {x} ^ {\prime} = 0, \quad x \in (0, 1), \quad t > 0, \\ u (x, 0) = 1 + \sin 2 \pi x, \quad x \in [ 0, 1 ], \\ u (1, t) = 1 + \sin 4 \pi t. \\ \end{array}
$$

此问题的精确解为 $u = 1 + \sin 2\pi (x + 2t)$

数值边值条件分别为

(a) $u_0^{n + 1} = u_0^n +\frac{2\tau}{h} (u_1^n -u_0^n);$   
(b) $u_0^n = u_1^n$   
(c) $u_0^{n + 1} - 2u_1^{n + 1} + u_2^{n + 1} = 0,$ 并与精确解比较

3. (实习题) 用差分格式

$$
\frac {u _ {j} ^ {n} - u _ {j} ^ {n - 1}}{\tau} + a \frac {u _ {j + 1} ^ {n - 1} - u _ {j - 1} ^ {n - 1}}{h} = \frac {u _ {j + 1} ^ {n - 1} - 2 u _ {j} ^ {n - 1} + u _ {j - 1} ^ {n - 1}}{h ^ {2}}
$$

和

$$
\frac {u _ {j} ^ {n} - u _ {j} ^ {n - 1}}{\tau} + a \frac {u _ {j} ^ {n - 1} - u _ {j - 1} ^ {n - 1}}{h} = \frac {u _ {j + 1} ^ {n - 1} - 2 u _ {j} ^ {n - 1} + u _ {j - 1} ^ {n - 1}}{h ^ {2}}
$$

求解

$$
\frac {\partial u}{\partial t} + b \frac {\partial u}{\partial x} = \frac {\partial^ {2} u}{\partial x ^ {2}}.
$$

初值条件：

$$
u (x, 0) = \left\{ \begin{array}{l l} x, & 0 \leqslant x \leqslant \frac {1}{2}, \\ (1 - x), & \frac {1}{2} \leqslant x \leqslant 1. \end{array} \right.
$$

边值条件：

$$
u (0, t) = u (1, t) = 0.
$$

取网比 $r = \frac{1}{2}$ , 试就 $b = 1,5,10,20$ 和步长 $h = 0.1,0.01$ , 计算 $t = 0.1,0.5$ 的差分解.

# 边值问题的变分形式与Ritz-Galerkin法

# 5.1 二次函数的极值

数学物理中的变分原理, 有重要的理论和实际意义, 也是构造微分方程数值解法的基础. 为了便于读者理解, 本节以二次函数的极值问题为例, 介绍变分原理的基本概念和方法. 主要结论和证明思想均可平行推广到线性对称微分方程边值问题.

在 $n$ 维欧氏空间 $\mathbb{R}^n$ 中引进向量、矩阵记号：

$$
\begin{array}{l} \boldsymbol {x} = \left(\xi_ {1}, \xi_ {2}, \dots , \xi_ {n}\right) ^ {\mathrm {T}}, \\ \boldsymbol {b} = \left(b _ {1}, b _ {2}, \dots , b _ {n}\right) ^ {\mathrm {T}}, \\ \boldsymbol {A} = \left[ \begin{array}{c c c c} a _ {1 1} & a _ {1 2} & \dots & a _ {1 n} \\ a _ {2 1} & a _ {2 2} & \dots & a _ {2 n} \\ \vdots & \vdots & & \vdots \\ a _ {n 1} & a _ {n 2} & \dots & a _ {n n} \end{array} \right]. \\ \end{array}
$$

这里以及往后仍用（） $\mathbf{\Psi}^{\mathrm{T}}$ 表示括号内向量或矩阵的转置. 令 $\pmb{y} = (\eta_{1}, \eta_{2}, \dots, \eta_{n})^{\mathrm{T}}$ , 定义 $\pmb{x}, \pmb{y}$ 的内积为

$$
\left(\boldsymbol {x}, \boldsymbol {y}\right) = \sum_ {i = 1} ^ {n} \xi_ {i} \eta_ {i}.
$$

考虑 $n$ 个变量的二次函数：

$$
F (\boldsymbol {x}) = F \left(\xi_ {1}, \xi_ {2}, \dots , \xi_ {n}\right) = \sum_ {i, j = 1} ^ {n} \alpha_ {i j} \xi_ {i} \xi_ {j} - \sum_ {i = 1} ^ {n} b _ {i} \xi_ {i} = (\boldsymbol {A x}, \boldsymbol {x}) - (\boldsymbol {b}, \boldsymbol {x}).
$$

它在 $\pmb{x}_0 = \left(\xi_1^{(0)},\xi_2^{(0)},\dots ,\xi_n^{(0)}\right)^{\mathrm{T}}$ 取极值的必要条件是

$$
\frac {\partial F \left(\xi_ {1} ^ {(0)} , \xi_ {2} ^ {(0)} , \cdots , \xi_ {n} ^ {(0)}\right)}{\partial \xi_ {k}} = \sum_ {i = 1} ^ {n} \left(\alpha_ {i k} + \alpha_ {k i}\right) \xi_ {i} ^ {(0)} - b _ {k} = 0, \quad k = 1, 2, \dots , n.
$$

假定 $\alpha_{ik} = \alpha_{ki}$ ，即 $\mathbf{A}$ 为对称矩阵，则

$$
2 \sum_ {i = 1} ^ {n} \alpha_ {k i} \xi_ {i} ^ {(0)} = b _ {k}, \quad k = 1, 2, \dots , n.
$$

不难看出，若令

$$
J (\boldsymbol {x}) = \frac {1}{2} (\boldsymbol {A x}, \boldsymbol {x}) - (\boldsymbol {b}, \boldsymbol {x}), \tag {5.1}
$$

则二次函数 $J(\pmb {x})$ 于 $\pmb{x}_0$ 取极值的必要条件是： $\pmb{x}_0$ 是线性代数方程组

$$
A x = b \tag {5.2}
$$

的解.

为了进一步研究 $J(\pmb {x})$ 于 $\pmb{x}_0$ 的极值性质，考虑实变量 $\lambda$ 的二次函数

$$
\varphi (\lambda) = J \left(\boldsymbol {x} _ {0} + \lambda \boldsymbol {x}\right),
$$

其中 $\pmb{x}$ 是任一 $n$ 维非零向量. 若 $J(x)$ 于 $\pmb{x}_0$ 取极小值, 则对任何 $\lambda \neq 0, \varphi(\lambda) = J(\pmb{x}_0 + \lambda \pmb{x}) > J(\pmb{x}_0) = \varphi(0)$ , 即 $\varphi(\lambda)$ 于 $\lambda = 0$ 取极小值. 反之, 若 $\varphi(\lambda)$ 于 $\lambda = 0$ 取极小值, 则对任何非零向量 $\pmb{x}, J(\pmb{x}_0 + \pmb{x}) = \varphi(1) > \varphi(0) = J(\pmb{x}_0)$ , 即 $J(\pmb{x})$ 于 $\pmb{x}_0$ 取极小值. 这样, 我们就把多变量函数 $J(\pmb{x})$ 的极值问题化成单变量函数 $\varphi(\lambda)$ 的极值问题

现在研究 $J(\pmb {x})$ 存在极小值的充要条件.显然

$$
\varphi (\lambda) = J \left(\boldsymbol {x} _ {0}\right) + \frac {\lambda}{2} \left[ \left(\boldsymbol {A} \boldsymbol {x} _ {0}, \boldsymbol {x}\right) + \left(\boldsymbol {A} \boldsymbol {x}, \boldsymbol {x} _ {0}\right) - 2 (\boldsymbol {b}, \boldsymbol {x}) \right] + \frac {\lambda^ {2}}{2} (\boldsymbol {A} \boldsymbol {x}, \boldsymbol {x}).
$$

因为 $A$ 是对称矩阵，故

$$
\varphi (\lambda) = J \left(\pmb {x} _ {0} + \lambda \pmb {x}\right) = J \left(\pmb {x} _ {0}\right) + \lambda \left(\pmb {A} \pmb {x} _ {0} - \pmb {b}, \pmb {x}\right) + \frac {\lambda^ {2}}{2} (\pmb {A} \pmb {x}, \pmb {x}). \tag {5.3}
$$

若 $J(\pmb {x})$ 于 $\pmb{x}_0$ 取极小值，则

$$
\varphi^ {\prime} (0) = \left(\boldsymbol {A} \boldsymbol {x} _ {0} - \boldsymbol {b}, \boldsymbol {x}\right) = 0, \quad \forall \boldsymbol {x} \in \mathbb {R} ^ {n},
$$

从而 $Ax_0 - b = 0$ ，这说明 $\pmb{x}_0$ 是(5.2)的解.又

$$
\varphi^ {\prime \prime} (0) = (\boldsymbol {A} \boldsymbol {x}, \boldsymbol {x}) > 0, \quad \forall \boldsymbol {x} \in \mathbb {R} ^ {n} \backslash \{\boldsymbol {0} \},
$$

故 $\mathbf{A}$ 必为正定矩阵

反之, 设 $\mathbf{A}$ 是对称正定矩阵, $\pmb{x}_0$ 是方程组 (5.2) 的解, 即

$$
A x _ {0} - b = 0,
$$

则由 (5.3) 得

$$
\varphi (\lambda) = J \left(\boldsymbol {x} _ {0}\right) + \frac {\lambda^ {2}}{2} (\boldsymbol {A x}, \boldsymbol {x}) = \varphi (0) + \frac {\lambda^ {2}}{2} (\boldsymbol {A x}, \boldsymbol {x}) > \varphi (0), \quad \lambda \neq 0, \quad \boldsymbol {x} \neq \mathbf {0}.
$$

这说明 $J(\pmb {x})$ 于 $\pmb{x}_0$ 取极小值.于是我们得

定理5.1 设矩阵 $\mathbf{A}$ 对称正定，则下列两问题等价：

(1) 求 $\pmb{x}_0 \in \mathbb{R}^n$ 使

$$
J \left(\boldsymbol {x} _ {0}\right) = \min  _ {\boldsymbol {x} \in \mathbb {R} ^ {n}} J (\boldsymbol {x}), \tag {5.4}
$$

其中 $J(\pmb {x})$ 是由(5.1)定义的二次函数

(2) 求下列方程组的解：

$$
\boldsymbol {A} \boldsymbol {x} = \boldsymbol {b}. \tag {5.5}
$$

$J(\pmb {x})$ 是定义在全空间 $\mathbb{R}^n$ 上的二次函数，称为 $\mathbb{R}^n$ 上的二次泛函或简称泛函数.泛函数 $J(\pmb {x})$ 由两部分组成：第一部分是二次项 $\frac{1}{2} (A\pmb {x},\pmb {x})$ ，它由矩阵 $\pmb{A}$ 决定；第二部分是一次项 $(b,x)$ ，它由右端向量 $\textit{\textbf{b}}$ 决定.

定理5.1表明，在矩阵 $\pmb{A}$ 为对称正定的条件下，若 $\pmb{x_0}$ 是极值问题(5.4）的解，则它也是线性方程组（5.5）的解；反之亦然.因此为了确定并计算 $\pmb{x_0}$ ，可采取两种不同途径：一种是求方程组（5.5）的解，另一种是求泛函数 $J(\pmb {x})$ 的极小值.我们更强调这后一途径因为许多数学物理问题，其直接的数学形式就是求意义更广的“二次泛函”的极小值，只是对解作了某些“光滑性”假设之后，才归结到微分方程.其次，即便是熟知的微分方程边值问题，我们也宁愿把它化为某一“二次泛函”的极小值问题，因为从极值问题出发建立数值解法往往更方便

# 5.1.1 习题

如果 $\varphi'(0) = 0$ ，则称 $\pmb{x}_0$ 是 $J(\pmb{x})$ 的驻点（或稳定点）。设矩阵 $\mathbf{A}$ 对称（不必正定），求证 $\pmb{x}_0$ 是 $J(\pmb{x})$ 的驻点的充要条件是： $\pmb{x}_0$ 是方程组 $\pmb{Ax} = \pmb{b}$ 的解。

# 5.2 Sobolev 空间初步

作为后面将要介绍的偏微分方程变分形式及多种数值方法的理论基础，本节我们引入 Sobolev 空间的概念并介绍一些基本性质（参见 [18, 24, 35]).

# 5.2.1 Sobolev空间的定义

由于偏微分方程涉及函数的导数，但经典导数的定义不能满足后面的理论需求，本小节我们先引入弱导数的概念，再结合我们学过的 $L^2$ 函数空间来引入 Sobolev 空间的定义.

弱导数

设 $\Omega \subset \mathbb{R}^d (d = 1,2,3)$ 是开集.记 $C_0^\infty (\varOmega)$ 是 $\varOmega$ 中的紧支集无穷次可微函数的集合定义 $L_{\mathrm{loc}}^{1}(\varOmega)$ 为局部可积函数的集合：

$$
L _ {\mathrm {l o c}} ^ {1} (\varOmega) = \left\{f: f | _ {K} \in L ^ {1} (K), \forall \text {紧 集} K \subset \varOmega \right\}.
$$

定义5.1（弱导数）假设 $f\in L_{\mathrm{loc}}^{1}(\Omega),1\leqslant i\leqslant d.$ 如果有 $g_{i}\in L_{\mathrm{loc}}^{1}(\Omega)$ 满足

$$
\int_ {\Omega} g _ {i} \varphi \mathrm {d} \boldsymbol {x} = - \int_ {\Omega} f \frac {\partial \varphi}{\partial x _ {i}} \mathrm {d} \boldsymbol {x}, \quad \forall \varphi \in C _ {0} ^ {\infty} (\Omega),
$$

那么称 $g_{i}$ 为 $f$ 关于 $x_{i}$ 在 $\Omega$ 上的弱（偏）导数，并记为

$$
\partial_ {x _ {i}} f = \frac {\partial f}{\partial x _ {i}} = g _ {i}, \quad i = 1, 2, \dots , d.
$$

类似地，对多重指标 $\alpha = (\alpha_{1},\alpha_{2},\dots ,\alpha_{d})\in \mathbb{N}^{d}$ ，记 $|\alpha | = \alpha_{1} + \alpha_{2} + \dots +\alpha_{d}$ ，可如下定义 $\partial^{\alpha}f\in L_{\mathrm{loc}}^{1}(\varOmega)$

$$
\int_ {\Omega} \partial^ {\alpha} f \varphi \mathrm {d} x = (- 1) ^ {| \alpha |} \int_ {\Omega} f \partial^ {\alpha} \varphi \mathrm {d} x, \quad \forall \varphi \in C _ {0} ^ {\infty} (\Omega),
$$

其中 $\partial^{\alpha} = \partial_{x_1}^{\alpha_1}\partial_{x_2}^{\alpha_2}\dots \partial_{x_d}^{\alpha_d}.$

显然弱导数概念是经典导数定义的推广，它保留了经典导数的分部积分的性质。但混合弱偏导数与求导次序无关。下面引理可以保证弱导数的唯一性。

引理5.1（变分学基本引理） 假设 $f\in L_{\mathrm{loc}}^{1}(\varOmega)$ 满足 $\int_{\varOmega}f\varphi \mathrm{d}\pmb {x} = 0,\forall \varphi \in C_0^\infty (\varOmega),$ 则 $f^{\mathrm{a.e.}} = 0$

例5.1 令 $d = 1, \Omega = (-1,1), f(x) = 1 - |x|$ . $f$ 的弱导数为

$$
g = \left\{ \begin{array}{l l} 1, & x \leqslant 0, \\ - 1, & x > 0. \end{array} \right.
$$

利用引理5.1可以证明 $g$ 的弱导数不存在，这里略去推导

Sobolev 空间

定义5.2(Sobolev空间） 对非负整数 $k$ ，定义

$$
H ^ {k} (\Omega) = \left\{u \in L ^ {2} (\Omega): \partial^ {\alpha} u \in L ^ {2} (\Omega), \quad \forall | \alpha | \leqslant k \right\}.
$$

可以证明 $H^{k}(\Omega)$ 在下面范数下是Banach空间：

$$
\| u \| _ {H ^ {k} (\Omega)} = \left(\sum_ {| \alpha | \leqslant k} \frac {| \alpha | !}{\alpha_ {1} ! \alpha_ {2} ! \cdots \alpha_ {n} !} \| \partial^ {\alpha} u \| _ {L ^ {2} (\Omega)} ^ {2}\right) ^ {\frac {1}{2}}.
$$

后面还会用到如下半范数的定义：

$$
| u | _ {H ^ {k} (\Omega)} = \left(\sum_ {| \alpha | = k} \frac {| \alpha | !}{\alpha_ {1} ! \alpha_ {2} ! \cdots \alpha_ {n} !} \| \partial^ {\alpha} u \| _ {L ^ {2} (\Omega)} ^ {p}\right) ^ {\frac {1}{2}}.
$$

后面往往简记 $\| \cdot \| _k = \| \cdot \|_{H^k (\Omega)},|\cdot |_k = |\cdot |_{H^k (\Omega)}.H^k (\Omega)$ 在下面内积下是Hilbert空间：

$$
(u, v) _ {k, \Omega} = \sum_ {| \alpha | \leqslant k} \frac {| \alpha | !}{\alpha_ {1} ! \alpha_ {2} ! \cdots \alpha_ {n} !} \int_ {\Omega} \partial^ {\alpha} u \partial^ {\alpha} v d x.
$$

记 $H_0^k (\varOmega)$ 为 $C_0^\infty (\varOmega)$ 在 $H^{k}(\varOmega)$ 中的闭包

例5.2 (1) 令 $\Omega = (0,1)$ , 考虑函数 $u = x^{\alpha}$ . 易知: 当 $\alpha > -\frac{1}{2}$ 时, $u \in L^{2}(\Omega)$ ; 当 $\alpha > \frac{1}{2}$ 或 $\alpha = 0$ 时, $u \in H^{1}(\Omega)$ ; 当 $\alpha > k - \frac{1}{2}$ 或 $\alpha = 0,1,\dots,k-1$ 时, $u \in H^{k}(\Omega)$ .

(2) 令 $\Omega = \left\{\pmb{x} \in \mathbb{R}^2 : |\pmb{x}| < \frac{1}{2}\right\}$ . 考虑函数 $f(\pmb{x}) = \ln |\ln |\pmb{x}|$ . 则 $f \in H^{1}(\Omega)$ , 但 $f \notin L^{\infty}(\Omega)$ . 此例说明, 在二维情形下, $H^{1}(\Omega)$ 中的函数可能不连续甚至无界.

# 5.2.2 Sobolev空间的性质

磨光

先介绍 Sobolev 空间中函数的磨光. 设 $\rho \in C_0^\infty(\mathbb{R}^d)$ 满足

$$
\rho (\boldsymbol {x}) \geqslant 0, \quad \int_ {\mathbb {R} ^ {d}} \rho (\boldsymbol {x})   d \boldsymbol {x} = 1, \quad \operatorname {s u p p} (\rho) \subset \{\boldsymbol {x}: | \boldsymbol {x} | \leqslant 1 \}.
$$

例如

$$
\rho (\boldsymbol {x}) = \left\{ \begin{array}{l l} C \mathrm {e} ^ {\frac {1}{| \boldsymbol {x} | ^ {2} - 1}}, & | \boldsymbol {x} | <   1, \\ 0, & | \boldsymbol {x} | \geqslant 1, \end{array} \right.
$$

其中常数 $C$ 使得 $\int_{\mathbb{R}^d}\rho (\pmb {x})\mathrm{d}\pmb {x} = 1.$ 对 $\varepsilon >0,$ 函数 $\rho_{\varepsilon}(\pmb {x}) = \varepsilon^{-d}\rho \left(\frac{\pmb{x}}{\varepsilon}\right)\in C_0^\infty (\mathbb{R}^d)$ 且 $\operatorname {supp}(\rho_{\varepsilon})\subset \{\pmb {x}:|\pmb {x}|\leqslant \varepsilon \} .$ 称 $\rho_{\varepsilon}$ 为磨光核 (mollifier）并称卷积

$$
u _ {\varepsilon} (\boldsymbol {x}) = \left(\rho_ {\varepsilon} * u\right) (\boldsymbol {x}) = \int_ {\mathbb {R} ^ {d}} \rho_ {\varepsilon} (\boldsymbol {x} - \boldsymbol {y}) u (\boldsymbol {y}) \mathrm {d} \boldsymbol {y}
$$

为 $u$ 的磨光函数（正则化，regularization）

引理5.2 (i) 如果 $u \in L_{\mathrm{loc}}^{1}(\mathbb{R}^{d})$ ，则 $\forall \varepsilon > 0, u_{\varepsilon} \in C^{\infty}(\mathbb{R}^{d})$ 且对任意多重指标 $\alpha$ 有 $\partial^{\alpha}(\rho_{\varepsilon} * u) = (\partial^{\alpha} \rho_{\varepsilon}) * u$ ；

(ii) 如果 $u \in L^{2}(\mathbb{R}^{d})$ , 则 $u_{\varepsilon} \in L^{2}(\mathbb{R}^{d}), \| u_{\varepsilon}\|_{L^{2}(\mathbb{R}^{d})} \leqslant \| u\|_{L^{2}(\mathbb{R}^{d})}$ , 且 $\lim_{\varepsilon \to 0}\| u_{\varepsilon} - u\|_{L^{2}(\mathbb{R}^{d})} = 0$

引理5.3 设 $\Omega$ 是一个区域（即连通的开集）， $u \in H^{1}(\Omega)$ ，且 $\nabla u \stackrel{\mathrm{a.e.}}{=} 0$ 于 $\Omega$ ，则 $u$ 在 $\Omega$ 上是常数。

证明 设 $K \coloneqq B(\pmb{x}_0, r) \subset \Omega$ 是任一球，取 $\varepsilon > 0$ 足够小使得 $B(\pmb{x}_0, r + \varepsilon) \subset \Omega$ . 将 $u$ 于 $\Omega$ 外作零延拓，仍记为 $u$ ，令 $u_{\varepsilon} \coloneqq \rho_{\varepsilon} * u$ . 显然， $\rho_{\varepsilon}(\pmb{x}, \cdot)|_{\Omega} \in C_0^\infty(\Omega)^d, \forall \pmb{x} \in K$ 故由弱偏导数的定义知，在 $K$ 中 $\nabla u_{\varepsilon} = (\nabla \rho_{\varepsilon}) * u = 0$ . 既然 $u_{\varepsilon}$ 是光滑的， $u_{\varepsilon}$ 在 $K$ 上是常数. 另外由引理5.2， $u_{\varepsilon} \to u$ 于 $L^2(K)$ . 故 $u$ 在 $K$ 上是常数. 由 $\Omega$ 的连通性及有限覆盖定理，得证. □

稠密性定理

定理5.2（稠密性定理）设 $k\geqslant 0,\Omega \subset \mathbb{R}^d$ 是有界多面体区域，则 $C^\infty (\bar{\Omega})\coloneqq$ $\{\varphi |_{\varOmega}:\varphi \in C_0^\infty (\mathbb{R}^d)\}$ 在 $H^{k}(\varOmega)$ 中稠密.

这里， $\mathbb{R}^d$ 中多面体在 $d = 1$ 时指线段， $d = 2$ 时指多边形， $d = 3$ 时就是通常的多面体。假设 $\Omega \subset \mathbb{R}^2$ 是单位圆，则存在线性算子 $\gamma: H^1(\Omega) \to L^2(\Omega)$ 。

嵌入定理

设 $X \subset Y$ 是两个Banach空间, 其范数分别记为 $\| \cdot \|_X, \| \cdot \|_Y$ . 我们称 $X$ 连续嵌入到 $Y$ , 记为 $X \hookrightarrow Y$ , 如果 $\| v \|_Y \leqslant C \| v \|_X, \forall v \in X$ . 我们称 $X$ 紧嵌入到 $Y$ , 记为 $X \hookrightarrow Y$ , 如果 $X$ 中的任意有界序列有在 $Y$ 中收敛的子序列. 显然当 $X \hookrightarrow Y$ 时, 恒等算子 $id_{X \leftrightarrow Y}$ 是连续的; 当 $X \hookrightarrow Y$ 时, 恒等算子 $id_{X \mapsto Y}$ 是紧算子; 并且 $X \hookrightarrow Y$ 蕴含 $X \hookrightarrow Y$ .

定理5.3（Sobolev嵌入定理） 假设 $\Omega \subset \mathbb{R}^d$ 是有界多面体区域， $k\geqslant 0$

(i) 如果 $k - \frac{d}{2} < 0$ ，则

$$
H ^ {k} (\varOmega) \hookrightarrow L ^ {q} (\varOmega), \quad \text {其 中} \quad q = {\frac {d}{d / 2 - k}},
$$

$$
H ^ {k} (\Omega) \hookrightarrow \hookrightarrow L ^ {q ^ {\prime}} (\Omega), \quad \forall 1 \leqslant q ^ {\prime} <   q.
$$

(ii) 如果 $k - \frac{d}{2} = 0$ ，则

$$
H ^ {k} (\Omega) \hookrightarrow \hookrightarrow L ^ {q} (\Omega), \quad \forall 1 \leqslant q <   \infty .
$$

(iii) 如果 $0 < \alpha \coloneqq k - \frac{d}{2} < 1$ , 则

$$
H ^ {k} (\Omega) \hookrightarrow C ^ {0, \alpha} (\bar {\Omega}), \quad H ^ {k} (\Omega) \hookrightarrow \hookrightarrow C ^ {0, \alpha^ {\prime}} (\bar {\Omega}), \quad \forall 0 \leqslant \alpha^ {\prime} <   \alpha .
$$

注5.1 若 $k - \frac{d}{2} \geqslant 1$ ，则可以对函数本身及导数应用上面的 Sobolev 嵌入定理，得到相应结果。

例5.3 设 $\Omega \subset \mathbb{R}^d$ 是有界多面体区域

(i) $H^{1}(\varOmega)\hookrightarrow\left\{\begin{array}{ll}C^{0,1/2}(\bar{\varOmega}),&d=1,\\ L^{q}(\varOmega),\quad 1\leqslant q< \infty,&d=2,\\ L^{6}(\varOmega),&d=3.\end{array}\right.$   
(ii) $H^{1}(\varOmega)\hookrightarrow\hookrightarrow L^{2}(\varOmega)$ ; 由归纳法得 $H^{k}(\varOmega)\hookrightarrow\hookrightarrow H^{k-1}(\varOmega), \forall k\geqslant 1$

定理5.4（Poincaré-Friedrichs不等式）假设 $\Omega \subset \mathbb{R}^d$ 是有界多面体区域，则

$$
\| u \| _ {L ^ {2} (\Omega)} \leqslant C \| \nabla u \| _ {L ^ {2} (\Omega)}, \quad \forall u \in H _ {0} ^ {1} (\Omega) \quad (\text {F r i e d r i c h s}),
$$

$$
\| u - u _ {\Omega} \| _ {L ^ {2} (\Omega)} \leqslant C \| \nabla u \| _ {L ^ {2} (\Omega)}, \quad \forall u \in H ^ {1} (\Omega) \quad \text {(P o i n c a r e)},
$$

其中 $u_{\Omega} = \frac{1}{|\Omega|}\int_{\Omega}u(x)\mathrm{d}x.$

证明 仅证第二个不等式. 记空间 $V = \{v \in H^{1}(\Omega) : v_{\Omega} = 0\}$ , 则 Poincaré 不等式等价于

$$
\| v \| _ {L ^ {2} (\Omega)} \leqslant C \| \nabla v \| _ {L ^ {2} (\Omega)}, \quad \forall v \in V.
$$

用反证法. 假设上面不等式不成立, 则存在序列 $\{v_{n}\} \subset V$ 使得

$$
\left\| v _ {n} \right\| _ {L ^ {2} (\Omega)} = 1, \quad \left\| \nabla v _ {n} \right\| _ {L ^ {2} (\Omega)} \leqslant \frac {1}{n}.
$$

由 $H^{1}(\Omega)\hookrightarrow \hookrightarrow L^{2}(\Omega)$ ，知存在子序列(仍记为） $v_{n}$ 及某函数 $v\in L^2 (\varOmega)$ 使得 $v_{n}\to v$ 于 $L^2 (\varOmega)$ .由弱导数的定义及在 $L^2 (\varOmega)^d$ 中 $\nabla v_{n}\rightarrow 0$ ，可得 $\nabla v = 0$ ，从而，由引理5.3知 $v$ 为常数,再由 $v_{\varOmega}=0$ 可得 $v = 0$ .这与 $\| v\|_{L^2 (\varOmega)} = 1$ 矛盾. □

迹定理

我们知道区域 $\Omega$ 上一个连续函数的边界值就是其在边界上的限制. 下面以 $H^{1}$ 空间为例, 讨论 Sobolev 空间中函数的“边界值”的定义. 由于 Sobolev 空间中两个仅在零测集上不同的函数被认为是同一个 (类) 函数, 所以对一个 Sobolev 空间中的函数, 允许任意改变其在某零测集上的值, 而被认为还是一个函数. 注意到相对于区域 $\Omega$ , 其边界 $\partial \Omega$ 正好是零测集, 那么该如何定义一个 Sobolev 空间中函数的边界值呢?

先考虑一个简单情形： $\Omega \subset \mathbb{R}^2$ 是单位圆盘，即

$$
\Omega = \left\{\boldsymbol {x}: | \boldsymbol {x} | <   1 \right\} = \{(r, \theta): 0 \leqslant r <   1, 0 \leqslant \theta <   2 \pi \}.
$$

设 $u\in C^1 (\bar{\Omega})$ ，考虑其在 $\partial \varOmega$ 上的限制：

$$
\begin{array}{l} u (1, \theta) ^ {2} = \int_ {0} ^ {1} \frac {\partial}{\partial r} \left(r ^ {2} u (r, \theta) ^ {2}\right) d r = \int_ {0} ^ {1} \left(2 r u ^ {2} + 2 r ^ {2} u u _ {r}\right) d r \\ \leqslant \int_ {0} ^ {1} \left(2 r u ^ {2} + 2 r ^ {2} | u | | \nabla u |\right) d r. \\ \end{array}
$$

从而

$$
\begin{array}{l} \| u \| _ {L ^ {2} (\partial \Omega)} ^ {2} = \int_ {0} ^ {2 \pi} u (1, \theta) ^ {2} \mathrm {d} \theta \leqslant \int_ {0} ^ {2 \pi} \int_ {0} ^ {1} 2 r \left(u ^ {2} + | u | | \nabla u |\right) \mathrm {d} r \mathrm {d} \theta \\ = 2 \int_ {\Omega} \left(u ^ {2} + | u | | \nabla u |\right) d x \leqslant 2 \| u \| _ {L ^ {2} (\Omega)} ^ {2} + 2 \| u \| _ {L ^ {2} (\Omega)} \| \nabla u \| _ {L ^ {2} (\Omega)} \\ = 2 \| u \| _ {L ^ {2} (\Omega)} \left(\| u \| _ {L ^ {2} (\Omega)} + \| \nabla u \| _ {L ^ {2} (\Omega)}\right) \\ \leqslant 2 \| u \| _ {L ^ {2} (\Omega)} \sqrt {2} \left(\| u \| _ {L ^ {2} (\Omega)} ^ {2} + \| \nabla u \| _ {L ^ {2} (\Omega)} ^ {2}\right) ^ {\frac {1}{2}} \\ = 2 \sqrt {2} \| u \| _ {L ^ {2} (\Omega)} \| u \| _ {H ^ {1} (\Omega)}. \\ \end{array}
$$

也就是, 当 $u \in C^{1}(\bar{\Omega})$ 时,

$$
\left\| u \right\| _ {L ^ {2} (\partial \Omega)} \leqslant \sqrt [ 4 ]{8} \left\| u \right\| _ {L ^ {2} (\Omega)} ^ {\frac {1}{2}} \left\| u \right\| _ {H ^ {1} (\Omega)} ^ {\frac {1}{2}}. \tag {5.6}
$$

下面证明 $H^{1}(\Omega)$ 中的函数可以定义边界值并满足上面不等式

引理5.4 假设 $\Omega \subset \mathbb{R}^2$ 是单位圆盘，则存在线性算子 $\gamma_0:H^1 (\varOmega)\mapsto L^2 (\partial \varOmega)$ 满足：

$$
\| \gamma_ {0} u \| _ {L ^ {2} (\partial \Omega)} \leqslant \sqrt [ 4 ]{8} \| u \| _ {L ^ {2} (\Omega)} ^ {\frac {1}{2}} \| u \| _ {H ^ {1} (\Omega)} ^ {\frac {1}{2}}, \quad \forall u \in H ^ {1} (\Omega).
$$

并且, 若 $u \in C^{1}(\bar{\Omega})$ , 则 $\gamma_{0} u = u|_{\partial \Omega}$ .

证明 由稠密性定理5.2, 存在 $u_{j} \in C^{1}(\bar{\Omega})$ 使得 $\| u - u_{j}\|_{H^{1}(\Omega)} \to 0, j \to \infty$ . 由(5.6)有

$$
\left\| u _ {j} - u _ {k} \right\| _ {L ^ {2} (\partial \Omega)} \leqslant \sqrt [ 4 ]{8} \left\| u _ {j} - u _ {k} \right\| _ {H ^ {1} (\Omega)}.
$$

所以, $u_{j}|_{\partial \Omega}$ 是 $L^2 (\partial \Omega)$ 中的Cauchy序列，从而收敛，定义其极限为 $\gamma_0 u$ ，即

$$
\gamma_ {0} u = \lim  _ {j \rightarrow \infty} u _ {j}.
$$

先证明此定义不依赖于序列 $\{u_j\}$ 的选取. 假设有另一个序列 $\{v_j\} \subset C^1 (\bar{\Omega})$ 满足 $\lim_{j\to \infty}\| u - v_j\|_{H^1 (\Omega)} = 0,$ 则由 (5.6)有

$$
\left\| u _ {j} - v _ {j} \right\| _ {L ^ {2} (\partial \Omega)} \leqslant \sqrt [ 4 ]{8} \left\| u _ {j} - v _ {j} \right\| _ {H ^ {1} (\Omega)} \rightarrow 0 \quad (j \rightarrow \infty).
$$

即 $\gamma_0 u$ 是唯一的. 于是有

$$
\begin{array}{l} \| \gamma_ {0} u \| _ {L ^ {2} (\partial \Omega)} = \lim  _ {j \rightarrow \infty} \| u _ {j} \| _ {L ^ {2} (\partial \Omega)} \leqslant \lim  _ {j \rightarrow \infty} \sqrt [ 4 ]{8} \| u _ {j} \| _ {L ^ {2} (\Omega)} ^ {\frac {1}{2}} \| u _ {j} \| _ {H ^ {1} (\Omega)} ^ {\frac {1}{2}} \\ = \sqrt [ 4 ]{8} \| u \| _ {L ^ {2} (\Omega)} ^ {\frac {1}{2}} \| u \| _ {H ^ {1} (\Omega)} ^ {\frac {1}{2}}. \\ \end{array}
$$

上面引理定义的 $\gamma_0 u$ 称为 $u$ 在 $\partial \Omega$ 上的“迹”，是光滑函数边界值的推广，也可称为 $u$ 的边界值。可以看出，一个 $H^1$ 空间中函数的边界值的定义采用了“稠密性论证”的流程，先取一列收敛于该函数的充分光滑函数序列，证明该序列的边界值序列收敛，定义其极限为该函数的边界值。更一般地，我们有下面的结论。

定理5.5 设 $\Omega \subset \mathbb{R}^d$ 是有界多面体区域，则存在有界线性算子 $\gamma_0: H^1(\Omega) \mapsto L^2(\partial \Omega)$ 及常数 $C > 0$ 使得

$$
\| \gamma_ {0} u \| _ {L ^ {2} (\partial \Omega)} \leqslant C \| u \| _ {L ^ {2} (\Omega)} ^ {\frac {1}{2}} \| u \| _ {H ^ {1} (\Omega)} ^ {\frac {1}{2}}, \quad \forall u \in H ^ {1} (\Omega).
$$

并且, 若 $u \in C^{1}(\bar{\Omega})$ , 则 $\gamma_{0} u = u|_{\partial \Omega}$ .

注5.2 (i)在定理5.5的条件下，我们有：

$$
H _ {0} ^ {1} (\Omega) = \left\{v \in H ^ {1} (\Omega): \gamma_ {0} v = 0 \right\}.
$$

(ii) $\gamma_0 u$ 可以理解为 $u$ 在 $\partial \Omega$ 上的值. 以后常常省略 $\gamma_0$ , 将 $\gamma_0 u$ 简记为 $u|_{\partial \Omega}$ .

定理5.6（Green第一公式）设 $\Omega \subset \mathbb{R}^d$ 是有界多面体区域， $\kappa \in (L^{\infty}(\varOmega))^{d\times d}$ ，V $\kappa \in (L^{\infty}(\varOmega))^{1\times d},u\in H^{2}(\varOmega),v\in H^{1}(\varOmega)$ ，则

$$
- \int_ {\Omega} \nabla \cdot (\kappa \nabla u) v d \boldsymbol {x} = \int_ {\Omega} \kappa \nabla u \cdot \nabla v d \boldsymbol {x} - \int_ {\partial \Omega} \kappa \nabla u \cdot \boldsymbol {n} v, \tag {5.7}
$$

其中 $\nabla \cdot \kappa$ 表示对 $\kappa$ 的列求散度所得行向量值函数, $n$ 是 $\partial \Omega$ 的单位外法向量, $\partial \Omega$ 上的积分在 $d = 2$ 和 $d = 3$ 时分别是第一型曲线和曲面积分.

证明 由稠密性论证, 只需对 $C^\infty(\bar{\Omega})$ 中的函数证明 Green 公式成立. 这可由恒等式

$$
\nabla \cdot (\kappa \nabla u) v + \kappa \nabla u \cdot \nabla v = \nabla \cdot ((\kappa \nabla u) v)
$$

及 Gauss 公式得到.

# 5.2.3 习题

1. 证明定理5.4中的Friedrichs不等式

# 5.3 两点边值问题

# 5.3.1 极小位能原理

考虑两点边值问题：

$$
L u = - \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) + q u = f, \quad x \in (a, b), \tag {5.8}
$$

$$
u (a) = 0, \quad u ^ {\prime} (b) = 0, \tag {5.9}
$$

其中 $p \in C^{1}(\bar{I})$ （一次连续可微函数空间）， $p(x) \geqslant \min_{x \in I} p(x) = p_{\min} > 0, q \in C(\bar{I}), q \geqslant 0, f \in H^{0}(I), \bar{I} = [a, b]$ ，构造泛函数

$$
\begin{array}{l} J (u) = \frac {1}{2} (L u, u) - (f, u) \\ = - \frac {1}{2} \int_ {a} ^ {b} \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) u \mathrm {d} x + \frac {1}{2} \int_ {a} ^ {b} q u ^ {2} \mathrm {d} x - \int_ {a} ^ {b} f u \mathrm {d} x. \\ \end{array}
$$

对右端第一项施行分部积分，并用边值条件(5.9)代入，得

$$
- \int_ {a} ^ {b} \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) u \mathrm {d} x = - p \frac {\mathrm {d} u}{\mathrm {d} x} u \Bigg | _ {a} ^ {b} + \int_ {a} ^ {b} p \left(\frac {\mathrm {d} u}{\mathrm {d} x}\right) ^ {2} \mathrm {d} x = \int_ {a} ^ {b} p \left(\frac {\mathrm {d} u}{\mathrm {d} x}\right) ^ {2} \mathrm {d} x.
$$

令

$$
a (u, v) = \int_ {a} ^ {b} \left(p \frac {\mathrm {d} u}{\mathrm {d} x} \frac {\mathrm {d} v}{\mathrm {d} x} + q u v\right) \mathrm {d} x, \tag {5.10}
$$

便得

$$
J (u) = \frac {1}{2} a (u, u) - (f, u). \tag {5.11}
$$

设 $H_E^1$ 为 $H^1$ 中满足左边值条件 $u(a) = 0$ 的函数组成的子空间. 考虑和 (5.8) (5.9) 相应的变分问题: 求 $u_* \in H_E^1$ , 使

$$
J \left(u _ {*}\right) = \min  _ {u \in H _ {E} ^ {1}} J (u). \tag {5.12}
$$

由(5.10)定义的 $a(u,v)$ 十分重要，它在今后的讨论中将起关键作用.显然 $a(u,v)$ 分别对 $u,v$ 都是线性泛函，即

$$
a \left(c _ {1} u _ {1} + c _ {2} u _ {2}, v\right) = c _ {1} a \left(u _ {1}, v\right) + c _ {2} a \left(u _ {2}, v\right),
$$

$$
a \left(u, c _ {1} v _ {1} + c _ {2} v _ {2}\right) = c _ {1} a \left(u, v _ {1}\right) + c _ {2} a \left(u, v _ {2}\right).
$$

因为 $c_{1}, c_{2}$ 是常数, 所以称为双线性泛函或双线性形式. 在讨论极值问题 (5.12) 之前, 我们先导出双线性形式 $a(u, v)$ 的几个基本性质.

首先 $a(u,v)$ 是对称形式，即

$$
a (u, v) = a (v, u), \quad {\text {对 任 意}}   u, v \in H ^ {1} (I).
$$

$a(u,v)$ 的对称性是由微分算子 $L$ 的对称性决定的.实际上，设 $u,v\in C^2 (I)$ ，且满足边值条件(5.9)，则

$$
\begin{array}{l} (L u, v) = \int_ {a} ^ {b} \left[ - \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) v + q u v \right] \mathrm {d} x \\ = \int_ {a} ^ {b} \left(p \frac {\mathrm {d} u}{\mathrm {d} x} \frac {\mathrm {d} v}{\mathrm {d} x} + q u v\right) \mathrm {d} x. \tag {5.13} \\ \end{array}
$$

对调 $u, v$ 后，等式右端不变，所以

$$
(L u, v) = (L v, u) = (u, L v). \tag {5.14}
$$

如此的 $L$ 称为对称算子. 其次

$$
a (u, u) = \int_ {a} ^ {b} \left[ p \left(\frac {\mathrm {d} u}{\mathrm {d} x}\right) ^ {2} + q u ^ {2} \right] \mathrm {d} x \geqslant p _ {\min } \int_ {a} ^ {b} \left(\frac {\mathrm {d} u}{\mathrm {d} x}\right) ^ {2} \mathrm {d} x. \tag {5.15}
$$

如果注意到任一 $u \in H_E^1$ 可表示为

$$
u (x) = \int_ {a} ^ {x} u ^ {\prime} (t) \mathrm {d} t,
$$

则由 Schwarz 不等式,

$$
\int_ {a} ^ {b} | u | ^ {2} \mathrm {d} x \leqslant \frac {1}{2} (b - a) ^ {2} \int_ {a} ^ {b} | u ^ {\prime} (t) | ^ {2} \mathrm {d} t,
$$

上面这种形式的不等式称为Poincaré-Friedrichs不等式(参见定理5.4).因而

$$
\begin{array}{l} \int_ {a} ^ {b} \left| u ^ {\prime} \right| ^ {2} \mathrm {d} x = \frac {1}{2} \left(\int_ {a} ^ {b} \left| u ^ {\prime} \right| ^ {2} \mathrm {d} x + \int_ {a} ^ {b} \left| u ^ {\prime} \right| ^ {2} \mathrm {d} x\right) \\ \geqslant \frac {1}{(b - a) ^ {2}} \int_ {a} ^ {b} | u | ^ {2} d x + \frac {1}{2} \int_ {a} ^ {b} | u ^ {\prime} | ^ {2} d x \\ \geqslant \bar {\gamma} \| u \| _ {1} ^ {2}, \tag {5.16} \\ \end{array}
$$

其中 $\bar{\gamma} = \min \left\{\frac{1}{2},\frac{1}{(b - a)^2}\right\} >0.$ 联立 (5.15)，(5.16) 并令 $\gamma = \bar{\gamma} p_{\mathrm{min}}$ ，得

$$
a (u, u) \geqslant \gamma \| u \| _ {1} ^ {2}, \quad \forall u \in H _ {E} ^ {1}. \tag {5.17}
$$

我们称满足不等式(5.17)的双线性形式为正定的或强制的.特别地，当 $u\in C^2 (\bar{I})$ 且满足边值条件（5.9）时，由（5.14）（5.17）得

$$
(L u, u) \geqslant \gamma \| u \| _ {1} ^ {2}.
$$

因此也说 $L$ 是正定算子

最后，由 Schwarz 不等式知 $a(u, v)$ 满足不等式

$$
\left| a (u, v) \right| \leqslant M \| u \| _ {1} \| v \| _ {1}, \quad u, v \in H ^ {1} (I), \tag {5.18}
$$

$M$ 是与 $u, v$ 无关的常数. 称 (5.18) 为连续性条件

现在回到变分问题 (5.12). 任取 $u \in H_E^1$ , 考虑实变量 $\lambda$ 的函数

$$
\begin{array}{l} \varphi (\lambda) = J \left(u _ {*} + \lambda v\right) \\ = \frac {1}{2} a (u _ {*} + \lambda v, u _ {*} + \lambda v) - (f, u _ {*} + \lambda v) \\ = \frac {1}{2} a (u _ {*}, u _ {*}) + \frac {\lambda}{2} [ a (u _ {*}, v) + a (v, u _ {*}) ] + \frac {\lambda^ {2}}{2} a (v, v) - (f, u _ {*}) - \lambda (f, v). \\ \end{array}
$$

由 $a(u_{*},v)$ 的对称性，得

$$
\varphi (\lambda) = J \left(u _ {*}\right) + \lambda \left[ a \left(u _ {*}, v\right) - (f, v) \right] + \frac {\lambda^ {2}}{2} a (v, v). \tag {5.19}
$$

今证下列变分原理

定理5.7 设 $f \in C(I), u_* \in C^2$ 是边值问题(5.8)(5.9)的解，则 $u_*$ 使 $J(u)$ 达到极小值. 反之，若 $u_* \in C^2 \cap H_E^1$ 使 $J(u)$ 达到极小值，则 $u_*$ 是边值问题(5.8)(5.9)的解

证明 注意当 $u_{*} \in C^{2} \cap H_{E}^{1}, v \in H_{E}^{1}$ 时，

$$
\begin{array}{l} a (u _ {*}, v) - (f, v) = \int_ {a} ^ {b} \left(p \frac {\mathrm {d} u _ {*}}{\mathrm {d} x} \frac {\mathrm {d} v}{\mathrm {d} x} + q u _ {*} v - f v\right) \mathrm {d} x \\ = p \frac {\mathrm {d} u _ {*}}{\mathrm {d} x} v \Bigg | _ {a} ^ {b} + \int_ {a} ^ {b} \left[ - \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u _ {*}}{\mathrm {d} x}\right) + q u _ {*} - f \right] v \mathrm {d} x \\ = \int_ {a} ^ {b} \left(L u _ {*} - f\right) v \mathrm {d} x + p (b) u _ {*} ^ {\prime} (b) v (b). \tag {5.20} \\ \end{array}
$$

如果 $u_{*}$ 是边值问题(5.8)(5.9)的解，则 $Lu_{*} - f = 0,u_{*}^{\prime}(b) = 0$ ，从而

$$
\varphi^ {\prime} (0) = a \left(u _ {*}, v\right) - (f, v) = 0, \quad \forall v \in H _ {E} ^ {1},
$$

注意 (5.20). 由 (5.19) 及 $a(u, v)$ 的正定性, 当 $\lambda \neq 0, v \neq 0$ 时,

$$
J \left(u _ {*} + \lambda v\right) = J \left(u _ {*}\right) + \frac {\lambda^ {2}}{2} a (v, v) > J \left(u _ {*}\right).
$$

这说明 $u_{*}$ 使 $J(u)$ 达到极小值

反之, 若 $u_*$ 使 $J(u)$ 达到极小值, 则由 (5.19) (5.20) 得

$$
\begin{array}{l} \varphi^ {\prime} (0) = a \left(u _ {*}, v\right) - (f, v) \\ = \int_ {a} ^ {b} (L u _ {*} - f) v d x + p (b) u _ {*} ^ {\prime} (b) v (b) = 0, \quad \forall v \in H _ {E} ^ {1}. \tag {5.21} \\ \end{array}
$$

特别地，取 $v\in C_0^\infty (I)$ ，则

$$
\int_ {a} ^ {b} \left(L u _ {*} - f\right) v d x = 0, \quad \forall v \in C _ {0} ^ {\infty} (I).
$$

根据变分法基本引理， $u_{*}$ 满足方程

$$
L u _ {*} - f = 0.
$$

于是 (5.21) 化为

$$
p (b) u ^ {\prime} _ {*} (b) v (b) = 0, \quad \forall v \in H _ {E} ^ {1}.
$$

注意 $p(b) > 0$ ，取 $v(x) = x - a,$ 则 $v\in H_E^1$ ，且 $v(b) > 0$ ，可见 $u_{*}$ 必须满足右边值条件

$$
u _ {*} ^ {\prime} (b) = 0.
$$

因为在力学、物理中, 二次泛函 $J(u)$ 表示能量, 所以也称定理5.7为极小位能原理. 应当注意的是, 我们仅就二次连续可微解 $u_{*}$ (称为古典解) 建立了边值问题和变分问题的等价性. 对于非光滑函数 $u_{*}$ , 说它是边值问题的解就没有意义了. 但是许多物理、力学现象, 必须用非光滑函数才能描述它. 比如前面所举弦平衡的例子, 若作用于弦的外力是集中荷载, 则弦的平衡曲线 $u = u(x)$ 不再有连续的二阶导数, 某些点甚至是没有导数的

“尖点”, 这时 $u(x)$ 在古典意义下不可能是 (5.8) (5.9) 的解. 而能量 $J(u)$ 的表达式是积分式 (5.15), 被积函数只含 $u$ 的一阶导数, 只要 $u$ 连续且按段连续可微, 则 $J(u)$ 有意义. 因此变分问题 (5.12) 允许非光滑解 $u_{*} = u_{*}(x)$ , 称之为两点边值问题 (5.8) (5.9) 的广义解或弱解.

边值问题可能有广义解但没有古典解. 定理 5.7 告诉我们, 当边值问题存在古典解时, 它一定是广义解. 反之, 若广义解存在且二次连续可微, 则广义解就是古典解.

按照变分法, 我们称 (5.8) 是和泛函数 $J(u)$ 相关的 Euler 方程

从定理5.7知道，左边值条件 $u(a) = 0$ 和右边值条件 $u^{\prime}(b) = 0$ 有重要区别。前者必须强加在变分问题所在的函数类上，称为强制边值条件或本质边值条件。后者不必对函数类作为条件提出，只要函数 $u_{*}(x)$ 使 $J(u)$ 取极小值，则它必然满足该条件，因此称为自然边值条件。在数值求解边值问题时，区别这两类条件很重要，这是从变分问题出发构造数值方法的一个优点。

![](images/2c9fccd524ca38f4e862dd97a9fe1cdb49395bee1596f4d8d964c3e9acd89f1f.jpg)

人物简介

# 5.3.2 虚功原理

以 $v$ 乘(5.8)两端，沿区间 $[a,b]$ 积分，得

$$
\int_ {a} ^ {b} (L u - f) v \mathrm {d} x = \int_ {a} ^ {b} \left[ - \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) v + q u v - f v \right] \mathrm {d} x = 0. \tag {5.22}
$$

利用分部积分和关于 $u, v$ 的边值条件 (5.9), 则

$$
\begin{array}{l} - \int_ {a} ^ {b} \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) v \mathrm {d} x = - p \frac {\mathrm {d} u}{\mathrm {d} x} v \Bigg | _ {a} ^ {b} + \int_ {a} ^ {b} p \frac {\mathrm {d} u}{\mathrm {d} x} \frac {\mathrm {d} v}{\mathrm {d} x} \mathrm {d} x \\ = \int_ {a} ^ {b} p \frac {\mathrm {d} u}{\mathrm {d} x} \frac {\mathrm {d} v}{\mathrm {d} x} \mathrm {d} x. \\ \end{array}
$$

以此代入到(5.22)式，得

$$
\int_ {a} ^ {b} \left(p \frac {\mathrm {d} u}{\mathrm {d} x} \frac {\mathrm {d} v}{\mathrm {d} x} + q u v - f v\right) \mathrm {d} x = 0.
$$

若注意到双线性形式 $a(u,v)$ 的表达式(5.10)，则上式可写成

$$
a (u, v) - (f, v) = 0. \tag {5.23}
$$

这也是边值问题 (5.8) (5.9) 的变分形式, 其确切提法将在定理 5.8 中给出.

对 $u\in C^2\cap H_E^1,v\in H_E^1$ ，根据(5.20)，方程(5.23)左端

$$
a (u, v) - (f, v) = \int_ {a} ^ {b} (L u - f) v d x + p (b) u ^ {\prime} (b) v (b).
$$

假若 $u$ 是边值问题(5.8)(5.9)的解，则对任意 $v\in H_E^1$ ， $u$ 满足(5.23).反之，若对任意

$v \in H_{E}^{1}, u \in H_{E}^{1}$ 满足 (5.23), 则可按定理 5.7 的证法, 推出 $u$ 是边值问题 (5.8) (5.9) 的解, 于是有

定理5.8 设 $u \in C^2$ ，则 $u$ 是边值问题（5.8）（5.9）的解的充要条件是： $u \in H_E^1$ 且满足变分方程

$$
a (u, v) - (f, v) = 0, \quad \forall v \in H _ {E} ^ {1}. \tag {5.24}
$$

在力学中, (5.24) 左端表示虚功, 所以也称定理 5.8 为虚功原理. 当 $u$ 是边值问题的古典解时, 它也是变分方程 (5.24) 的解. 像位能原理一样, 变分方程 (5.24) 还允许非古典解, 我们称这样的解为边值问题的广义解或弱解.

虚功原理比位能原理更具有一般性，它不仅适用于对称正定算子方程（相当于力学中的保守场方程），而且也适用于非对称正定算子方程（非保守场方程）。实际上，定理5.8可直接推广为

定理5.9 设 $u \in C^2$ ，则 $u$ 满足

$$
\left\{ \begin{array}{l} L u = - \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) + r \frac {\mathrm {d} u}{\mathrm {d} x} + q u = f, \\ u (a) = 0, \quad u ^ {\prime} (b) = 0 \end{array} \right.
$$

的充要条件是 $u \in C^2 \cap H_E^1$ 且满足变分方程：

$$
\left\{ \begin{array}{l} a (u, v) - (f, v) = 0, \quad \forall v \in H _ {E} ^ {1}, \\ a (u, v) = \int_ {a} ^ {b} \left(p \frac {\mathrm {d} u}{\mathrm {d} x} \frac {\mathrm {d} v}{\mathrm {d} x} + r \frac {\mathrm {d} u}{\mathrm {d} x} v + q u v\right) \mathrm {d} x, \end{array} \right.
$$

其中 $p\in C^1,p_{\min} > 0,r,q\in C,f\in L^2$

此时双线性形式 $a(u,v)$ 非对称正定，除非 $r\equiv 0,q\geqslant 0$

# 5.3.3 习题

1. 证明非齐次两点边值问题

$$
\left\{ \begin{array}{l} L u = - \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) + q u = f, \quad a <   x <   b, \\ u (a) = \alpha , \quad u ^ {\prime} (b) = \beta \end{array} \right. \tag {5.25}
$$

与下列变分问题等价：求 $u_{*}\in H^{1},u_{*}(a) = \alpha$ ，使

$$
J(u_{*}) = \min_{\substack{u\in H^{1}\\ u(a) = \alpha}}J(u),
$$

其中

$$
J (u) = \frac {1}{2} a (u, u) - (f, u) - p (b) \beta u (b),
$$

而 $a(u,v)$ 如（5.10）（提示：先把边值条件齐次化）

2. 就边值问题 (5.25) 建立虚功原理  
3. 试建立与边值问题

$$
\left\{ \begin{array}{l l} L u = \frac {\mathrm {d} ^ {4} u}{\mathrm {d} x ^ {4}} + u = f, & a <   x <   b, \\ u (a) = u ^ {\prime} (a) = 0, & u (b) = u ^ {\prime} (b) = 0 \end{array} \right.
$$

等价的变分问题

# 5.4 二阶椭圆边值问题

为了记号简单, 我们用黑体的坐标向量 $\pmb{x} \in \mathbb{R}^2$ 表示 $\pmb{x} = (x,y)$ .

# 5.4.1 极小位能原理

作为模型，考虑 Poisson 方程的第一边值问题：

$$
- \Delta u = f (\boldsymbol {x}), \quad \boldsymbol {x} \in \Omega , \tag {5.26}
$$

$$
u | _ {\Gamma} = 0, \tag {5.27}
$$

其中边界 $\Gamma$ 为分段光滑曲线.作泛函数

$$
\begin{array}{l} J (u) = \frac {1}{2} (- \Delta u, u) - (f, u) \\ = \frac {1}{2} \int_ {\Omega} (- \Delta u) u \mathrm {d} x - \int_ {\Omega} f u \mathrm {d} x. \tag {5.28} \\ \end{array}
$$

利用Green第一公式(5.7)，我们得

$$
\int_ {\Omega} (- \Delta u) v \mathrm {d} x = \int_ {\Omega} \nabla u \cdot \nabla v \mathrm {d} x - \int_ {\Gamma} \frac {\partial u}{\partial n} v \mathrm {d} s, \tag {5.29}
$$

其中 $\pmb{n}$ 表示曲线边界 $\Gamma$ 的单位外法向量， $\frac{\partial u}{\partial \pmb{n}}$ 是 $u$ 沿 $\pmb{n}$ 的方向导数. 若 $u, v$ 满足边值条件 (5.27)，则

$$
\int_ {\Omega} (- \Delta u) v \mathrm {d} x = \int_ {\Omega} \nabla u \cdot \nabla v \mathrm {d} x. \tag {5.30}
$$

定义双线性形式

$$
a (u, v) = \int_ {\Omega} \nabla u \cdot \nabla v \mathrm {d} x. \tag {5.31}
$$

由(5.28)(5.30)和(5.31)，可将泛函数 $J(u)$ 写成

$$
J (u) = \frac {1}{2} a (u, u) - (f, u). \tag {5.32}
$$

在力学上， $J(u)$ 表示位能

从 (5.31), (5.32) 知道, 只要 $u \in H^{1}(\Omega), f \in L^{2}(\Omega)$ , 则 $J(u)$ 有意义. 此外还要求 $u$ 满足第一边值条件 (5.27). 以下用 $H_{0}^{1}(\Omega)$ 表示 $H^{1}(\Omega)$ 中一切满足 (5.27) 的函数组成的子空间.

现在提如下变分问题：求 $u_{*}\in H_{0}^{1}(\varOmega)$ ，使

$$
J (u _ {*}) = \min  _ {u \in H _ {0} ^ {1} (\Omega)} J (u). \tag {5.33}
$$

为了建立边值问题（5.26）（5.27）和变分问题（5.33）的等价性，先讨论双线性形式 $a(u,v)$ 的两个基本性质.

(1) 对称性. 显然

$$
a (u, v) = a (v, u), \quad \forall u, v \in H ^ {1} (\Omega).
$$

(2) 正定性. 对于 $u \in H_0^1(\Omega)$ ,

$$
a (u, u) = \int_ {\Omega} | \nabla u | ^ {2} \mathrm {d} x \mathrm {d} y = \| \nabla u \| ^ {2}.
$$

由Poincaré-Friedrichs不等式（见定理5.4）得：

$$
\| \nabla u \| ^ {2} \geqslant \gamma \| u \| _ {1} ^ {2}, \quad u \in H _ {0} ^ {1} (\Omega), \tag {5.34}
$$

其中 $\gamma > 0$ 是和 $u$ 无关的常数. 于是

$$
a (u, u) \geqslant \gamma \| u \| _ {1} ^ {2}, \quad u \in H _ {0} ^ {1} (\Omega). \tag {5.35}
$$

这说明 $a(u,v)$ 正定

其次, 对于 $u, v \in H^{1}(\Omega), a(u,v)$ 满足不等式 (连续性条件)

$$
| a (u, v) | \leqslant \| u \| _ {1} \| v \| _ {1}. \tag {5.36}
$$

由于 $a(u,v)$ 对称正定，也称 $-\Delta$ 为对称正定算子（参看(5.30)(5.31)).

对于 $u_{*}, u \in H_{0}^{1}(\Omega)$ , 考虑实参数 $\lambda$ 的函数

$$
\varphi (\lambda) = J \left(u _ {*} + \lambda u\right).
$$

利用 $a(u,v)$ 的对称性，可知

$$
\varphi (\lambda) = J \left(u _ {*}\right) + \lambda \left[ a \left(u _ {*}, u\right) - (f, u) \right] + \frac {\lambda^ {2}}{2} a (u, u). \tag {5.37}
$$

它和 (5.19) 有完全相同的形式

若进一步假定 $u_{*}\in C^{2}(\bar{\Omega})\cap H_{0}^{1}(\Omega)$ ，则由（5.30）（5.31）得出

$$
a \left(u _ {*}, u\right) - (f, u) = (- \Delta u _ {*} - f, u). \tag {5.38}
$$

设 $u_{*}$ 是边值问题(5.26)(5.27）的解，则

$$
\varphi^ {\prime} (0) = a \left(u _ {*}, u\right) - (f, u) = (- \Delta u _ {*} - f, u) = 0, \quad \forall u \in H _ {0} ^ {1} (\Omega),
$$

从而

$$
J \left(u _ {*} + \lambda u\right) = J \left(u _ {*}\right) + \frac {\lambda^ {2}}{2} a (u, u) > J \left(u _ {*}\right),
$$

对任意 $u \neq 0, u \in H_0^1(\Omega), \lambda \neq 0$ . 这说明 $u_*$ 使 $J(u)$ 达到极小.

与定理5.7类似，可证明使 $J(u)$ 达到极小的 $u_{*}$ ，当其属于 $C^2 (\bar{\Omega})\cap H_0^1 (\Omega)$ 时，必为(5.26）（5.27）的解.于是得

定理5.10 设 $u_* \in C^2(\bar{\Omega})$ 是边值问题(5.26)(5.27)的解，则 $u_*$ 使 $J(u)$ 达到极小。反之，若 $u_* \in C^2(\bar{\Omega}) \cap H_0^1(\Omega)$ 使 $J(u)$ 达到极小，则 $u_*$ 是边值问题(5.26)(5.27)的解（古典解）。

由于 $J(u)$ 在力学、物理学中表示能量，所以也称定理5.10为极小位能原理。注意定理5.10要求 $u_{*} \in C^{2}(\bar{\Omega})$ ，而变分问题（5.33）还允许不属于 $C^{2}(\bar{\Omega})$ 的解，称之为边值问题的广义解。

注5.3 若代替（5.27）的是非齐边值条件

$$
\left. u \right| _ {\Gamma} = \varphi (\boldsymbol {x}), \quad \varphi \in C ^ {1} (\Gamma), \tag {5.39}
$$

则取一特定函数 $u_0 \in C^2(\bar{\Omega})$ 满足 $u_0|_{\Gamma} = \varphi$ . 令 $v = u - u_0$ , 则 $v$ 满足方程:

$$
- \Delta v = f + \Delta u _ {0}
$$

和齐次边值条件 (5.27). 构造 $v$ 的二次泛函

$$
\begin{array}{l} \widehat {J} (v) = \frac {1}{2} \int_ {\Omega} | \nabla v | ^ {2} d x - \int_ {\Omega} \left(f + \Delta u _ {0}\right) v d x \\ = \frac {1}{2} \int_ {\Omega} | \nabla (u - u _ {0}) | ^ {2} d x - \int_ {\Omega} (f + \Delta u _ {0}) (u - u _ {0}) d x \\ = \frac {1}{2} \int_ {\Omega} \left(| \nabla u | ^ {2} - 2 f u\right) \mathrm {d} x - \int_ {\Omega} \nabla u \cdot \nabla u _ {0} \mathrm {d} x - \int_ {\Omega} \Delta u _ {0} u \mathrm {d} x + \text {常 数}. \\ \end{array}
$$

由 (5.29)

- $\int_{\Omega}\Delta u_0u\mathrm{d}\pmb {x} - \int_{\Omega}\nabla u\cdot \nabla u_0\mathrm{d}\pmb {x} = -\int_{\Gamma}u\frac{\partial u_0}{\partial\pmb{n}}\mathrm{d}s = -\int_{\Gamma}\varphi \frac{\partial u_0}{\partial\pmb{n}}\mathrm{d}s =$ 常数，足见

$$
\widehat {J} (v) = J (u) + \text {常 数}.
$$

由此可见，变分问题

$$
\widehat{J} (v_{*}) = \min_{v\in H^{1}_{0}}\widehat{J} (v)
$$

和

$$
J (u _ {*}) = \min  _ { \begin{array}{c} u \in H ^ {1} (\Omega) \\ u | _ {\Gamma} = \varphi \end{array} } J (u) \tag {5.40}
$$

等价, 且 $v_{*} = u_{*} - u_{0}$

据定理5.10，非齐次边值问题(5.26)(5.39)与变分问题(5.40）等价

# 5.4.2 自然边值条件

考虑第二、第三边值条件

$$
\left. \frac {\partial u}{\partial \boldsymbol {n}} + \alpha u \right| _ {\Gamma} = 0, \quad \alpha \geqslant 0. \tag {5.41}
$$

利用公式 (5.29),

$$
\begin{array}{l} J (u) = \frac {1}{2} (- \Delta u, u) - (f, u) \\ = \frac {1}{2} \int_ {\Omega} | \nabla u | ^ {2} \mathrm {d} x - \frac {1}{2} \int_ {\Gamma} \frac {\partial u}{\partial \boldsymbol {n}} u \mathrm {d} s - \int_ {\Omega} f u \mathrm {d} x \\ = \frac {1}{2} \int_ {\Omega} | \nabla u | ^ {2} \mathrm {d} x + \frac {1}{2} \int_ {\Gamma} \alpha u ^ {2} \mathrm {d} s - \int_ {\Omega} f u \mathrm {d} x. \\ \end{array}
$$

令

$$
a (u, v) = \int_ {\Omega} \nabla u \cdot \nabla v d x + \int_ {\Gamma} \alpha u v d s, \tag {5.42}
$$

则

$$
J (u) = \frac {1}{2} a (u, u) - (f, u). \tag {5.43}
$$

设 $u_{*}\in C^{2}(\bar{\varOmega}),u\in H^{1}(\varOmega)$ 由公式(5.29)，我们有

$$
a \left(u _ {*}, u\right) - (f, u) = (- \Delta u _ {*} - f, u) + \int_ {\Gamma} \left(\alpha u _ {*} + \frac {\partial u _ {*}}{\partial \boldsymbol {n}}\right) u d s. \tag {5.44}
$$

今考虑实变量 $\lambda$ 的函数

$$
\varphi (\lambda) = J \left(u _ {*} + \lambda u\right), \quad u _ {*}, u \in H ^ {1}.
$$

直接计算，可得形如（5.37）的展开式：

$$
\varphi (\lambda) = J \left(u _ {*}\right) + \lambda \left[ a \left(u _ {*}, u\right) - (f, u) \right] + \frac {\lambda^ {2}}{2} a (u, u), \tag {5.45}
$$

其中 $a(u_{*},u)$ 由(5.42)定义.与定理5.7的证法类似，可以证明

定理5.11 边值问题（5.26）（5.41）的解 $u_{*} \in C^{2}(\bar{\Omega})$ 是下列变分问题的解：求 $u_{*} \in H_{0}^{1}(\Omega)$ ，使

$$
J \left(u _ {*}\right) = \min  _ {u \in H ^ {1}} J (u). \tag {5.46}
$$

反之，变分问题(5.46）的解 $u_{*}$ 若属于 $C^2 (\bar{\Omega})$ ，则也是边值问题(5.26)(5.41）的解

若 $u_{*}\in H^{1}(\bar{\Omega})$ 是(5.46)的解，则称之为边值问题的广义解

值得指出的是，变分问题(5.46)并不要求 $u$ 满足任何边值条件，而它的解 $u_{*}$ 却自动满足(5.41)，这是第二、三边值条件与第一边值条件的一个重大差别.像两点边值问题一样，我们称第一边值条件为本质边值条件，第二、三边值条件为自然边值条件

# 5.4.3 虚功原理

像本章第3节那样，同样可以建立第一、第二、第三边值问题的虚功原理.为叙述统一，我们考虑Poisson方程(5.26）的混合边值问题.如图5.1，设边界 $\varGamma$ 分成互不相交的两部分： $\varGamma_{1}$ 和 $\varGamma_{2}$ .在 $\varGamma_{1}$ 上满足第一边值条件：

$$
u \mid_ {\Gamma_ {1}} = 0, \quad (u, u \Delta -) \frac {1}{c} = (u) L \tag {5.47}
$$

在 $\Gamma_{2}$ 上满足第二或第三边值条件：

$$
\left. \frac {\partial u}{\partial \boldsymbol {n}} + \alpha u \right| _ {\Gamma_ {2}} = 0, \quad \alpha \geqslant 0. \tag {5.48}
$$

![](images/59a414bfb058e9b325a7a2c004a6f28002aaf54095679d6b6fb3024852ce31c6.jpg)  
图5.1 边界

以 $H_E^1 (\Omega)$ 表示 $H^{1}(\Omega)$ 中满足第一边值条件(5.47）的函数组成的子空间．以 $\boldsymbol {v}\in$ $H_E^1 (\Omega)$ 乘（5.26）两端并在 $\varOmega$ 上积分，得

$$
\int_ {\Omega} [ (- \Delta u) v - f v ] \mathrm {d} x = 0. \tag {5.49}
$$

利用公式 (5.29) 及关于 $u, v$ 的边值条件 (5.47) (5.48) 得

$$
\begin{array}{l} \int_ {\Omega} (- \Delta u) v \mathrm {d} x \mathrm {d} y = \int_ {\Omega} \nabla u \cdot \nabla v \mathrm {d} x - \int_ {\Gamma} \frac {\partial u}{\partial n} v \mathrm {d} s \\ = \int_ {\Omega} \nabla u \cdot \nabla v \mathrm {d} x + \int_ {\Gamma_ {2}} \alpha u v \mathrm {d} s. \tag {5.50} \\ \end{array}
$$

定义双线性形式

$$
a (u, v) = \int_ {\Omega} \nabla u \cdot \nabla v d x + \int_ {\Gamma_ {2}} a u v d s, \tag {5.51}
$$

则可将(5.49)写成

$$
a (u, v) - (f, v) = 0. \tag {5.52}
$$

今提如下变分问题: 求 $u \in H_E^1(\Omega)$ , 使 $u$ 对一切 $v \in H_E^1(\Omega)$ 满足 (5.52).

设 $u\in C^2 (\bar{\Omega}),v\in H_E^1 (\Omega)$ ，则由(5.29）得

$$
a (u, v) - (f, v) = \int_ {\Omega} (- \Delta u - f) v d x + \int_ {\Gamma_ {2}} \left(\frac {\partial u}{\partial n} + \alpha u\right) v d s. \tag {5.53}
$$

与定理5.7的证法类似，可推出

定理5.12 设 $u \in C^2(\bar{\Omega})$ ，则 $u$ 满足(5.26)(5.47)(5.48)的充要条件是： $u \in H_E^1$ 且对任意 $v \in H_E^1$ 满足变分方程(5.52).

因为(5.52)左端在力学中表示虚功，故亦称定理5.12为虚功原理。和边值问题不同，变分方程（5.52）允许有不属于 $C^2(\Omega)$ 的解，称为边值问题的广义解或弱解。从定理5.12看出，边值条件（5.47）和（5.48）有重要差别，前者为本质边值条件，后者为自然边值条件。正如本章第3节指出过的，虚功原理较极小位能原理应用更广，它不必要求边值问题对称正定。

# 5.4.4 习题

1. 设 $u \in C(\bar{\Omega})$ 满足

$$
\int_ {\Omega} u \varphi \mathrm {d} x = 0, \quad \forall \varphi \in C _ {0} ^ {\infty} (\Omega),
$$

试证 $u\equiv 0$

2. 证明定理5.10的第二部分  
3. 试就 Poisson 方程 (5.26) 的非齐次边值条件

$$
\frac {\partial u}{\partial \boldsymbol {n}} + \alpha u | _ {\Gamma} = \beta , \quad \alpha \geqslant 0, \tag {5.54}
$$

导出等价的变分问题

4. 试就椭圆型方程第一边值问题：

$$
- \nabla \cdot (\kappa \nabla u) + \sigma u = f, \quad (x, y) \in \Omega , \quad u | _ {\Gamma} = g \tag {5.55}
$$

建立等价的极小位能原理和虚功原理, 其中 $\kappa = \kappa (x,y)\in C^{1}(\bar{\Omega}),\min_{\bar{\Omega}}\kappa >0,\sigma \in$ $C(\bar{\Omega}),\sigma \geqslant 0,f\in L^{2}(\varOmega),g\in C^{1}(\varGamma)$ ，而

$$
\nabla \cdot (\kappa \nabla u) = \frac {\partial}{\partial x} \left(\kappa \frac {\partial u}{\partial x}\right) + \frac {\partial}{\partial y} \left(\kappa \frac {\partial u}{\partial y}\right).
$$

# 5.5 Ritz-Galerkin 方法

前面各节讨论了如何化边值问题为等价的变分问题，本节讨论如何解相应的变分问题.必须指出，除少数特殊情形外，一般不可能求得问题的准确解，因此需要各种近似或数值解法.Ritz-Galerkin方法是最重要的一种解法，它是以后要讨论的有限元法的基础

用 $U$ 表示 $H_0^1, H_E^1, H^1$ 等 Sobolev 空间, $H = H^0$ 是 $L^2$ 空间. $L$ 代表本章第 3 节和第 4 节中的微分算子 (二阶常微分或偏微分算子). $a(u,v)$ 是由 $L$ 及边值条件决定的双线性形式, 它由 $(Lu,v)$ 经过分部积分并代入边值条件后得到. 得出 $a(u,v)$ 的表达式后, $u, v$ 就无需满足自然边值条件了, 但本质边值条件仍需满足, 就是说, $u, v$ 应属于空间 $U$ . 前已证明, $a(u,v)$ 是对称正定双线性形式, 即满足

$$
a (u, v) = a (v, u), \quad \forall u, v \in U, \tag {5.56}
$$

$$
a (u, u) \geqslant \gamma \| u \| _ {1} ^ {2}, \quad \forall u \in U, \tag {5.57}
$$

其中 $\gamma > 0$ 是与 $u$ 无关的常数. 正定性 (5.57) 也通常被称为强制性.

此外， $a(u,v)$ 还满足连续性条件

$$
\left| a (u, v) \right| \leqslant M \| u \| _ {1} \| v \| _ {1}, \quad u, v \in U, \tag {5.58}
$$

参看 (5.18) (5.36).

设 $f\in H$ ，则本章第1、3和4节的二次泛函可统一写成形式：

$$
J (u) = \frac {1}{2} a (u, u) - (f, u).
$$

于是边值问题 $Lu = f$ 等价于求 $u \in U$ , 使

$$
J (u) = \min  _ {v \in U} J (v). \tag {5.59}
$$

这就是极小位能原理

边值问题的另一变分形式是：求 $u\in U$ ，使

$$
a (u, v) = (f, v), \quad \forall v \in U. \tag {5.60}
$$

这就是虚功原理. 虚功原理并不要求 $a(u, v)$ 对称正定

变分问题 (5.59) 或 (5.60) 的主要困难是在无穷维空间 $U$ 上求解. Ritz-Galerkin 方法的基本思想在于用有穷维空间近似代替无穷维空间, 从而化成在有限维空间上近似求解 (参看本章第 1 节). 关键是如何选取有穷维空间.

设 $U_{n}$ 是 $U$ 的 $n$ 维子空间. 将极小位能原理 (5.59) 中的无穷维空间 $U$ 换为 $U_{n}$ 就得到求解边值问题的Ritz方法: 求 $u_{n} \in U_{n}$ 使得

$$
J \left(u _ {n}\right) = \min  _ {v _ {n} \in U _ {n}} J \left(v _ {n}\right). \tag {5.61}
$$

既然是在 $U_{n}$ 中找近似解， $U_{n}$ 也称为试探函数空间

类似地，如果将虚功原理即变分问题（5.60）中 $U$ 替换为 $U_{n}$ ，那么就得到求解边值问题的 Galerkin 方法：求 $u_{n} \in U_{n}$ 使得

$$
a \left(u _ {n}, v _ {n}\right) = (f, v _ {n}), \quad \forall v _ {n} \in U _ {n}. \tag {5.62}
$$

为了程序实现, 我们给出两种方法的矩阵形式. 设 $\varphi_1, \varphi_2, \dots, \varphi_n$ 是 $U_n$ 的一组基底, 称为基函数, 则 $U_n$ 中任一函数 $u_n$ 可表示为

$$
u _ {n} = \sum_ {j = 1} ^ {n} c _ {j} \varphi_ {j}. \tag {5.63}
$$

注意

$$
J (u _ {n}) = \frac {1}{2} a (u _ {n}, u _ {n}) - (f, u _ {n}) = \frac {1}{2} \sum_ {i, j = 1} ^ {n} a (\varphi_ {i}, \varphi_ {j}) c _ {i} c _ {j} - \sum_ {j = 1} ^ {n} c _ {j} (f, \varphi_ {j})
$$

是 $c_{1}, c_{2}, \dots, c_{n}$ 的二次函数， $a(\varphi_i, \varphi_j) = a(\varphi_j, \varphi_i)$ . 记

$$
\boldsymbol {A} = \left(a \left(\varphi_ {j}, \varphi_ {i}\right)\right) _ {n \times n}, \quad \boldsymbol {c} = \left(c _ {i}\right) _ {n \times 1}, \quad \boldsymbol {b} = \left(\left(f, \varphi_ {i}\right)\right) _ {n \times 1}, \quad I (\boldsymbol {c}) = \frac {1}{2} (\boldsymbol {A c}, \boldsymbol {c}) - (\boldsymbol {b}, \boldsymbol {c}).
$$

则显然 $J(u_{n}) = I(c)$ .Ritz法可改写为：求 $\pmb {c}\in \mathbb{R}^n$ 使得

$$
I (\boldsymbol {c}) = \min  _ {\boldsymbol {x} \in \mathbb {R} ^ {n}} I (\boldsymbol {x}). \tag {5.64}
$$

类似地, 在 (5.62) 中取 $v = \varphi_{i}, 1 \leqslant i \leqslant n$ , 知 Galerkin 法可改写为: 求 $c \in \mathbb{R}^{n}$ 使得

$$
\boldsymbol {A} \boldsymbol {c} = \boldsymbol {b}. \tag {5.65}
$$

显然解此方程组得到 $c$ 代入 (5.63) 就可以得到 Galerkin 法的解 $u_{n}$ . 注意到 $(A c, c) = a(u_{n}, u_{n})$ , 易知, 如果双线性形式 $a(u, v)$ 对称正定, 即满足 (5.56) (5.57), 那么 $A$ 是对称正定矩阵. 由此即知 (5.65) 唯一可解. 进一步由定理 5.1 知 (5.64) 与 (5.65) 等价, 故 Ritz 方法 (5.61) 和 Galerkin 法 (5.62) 也等价. 总结一下, 我们有如下定理:

定理5.13 假设双线性形式 $a(u,v)$ 对称正定，即满足(5.56)(5.57)，那么Ritz方法（5.61）和Galerkin法（5.62）等价且唯一可解

由此等价性, 习惯上称 (5.65) 为 Ritz-Galerkin 方程. 尽管 Ritz 法和 Galerkin 法导出的近似解 $u_{n}$ 及计算方法完全一样, 但二者的基础不同. Ritz 法基于极小位能原理, 而 Galerkin 法基于虚功原理, 所以 Galerkin 法较 Ritz 法应用更广, 方法推导也更直接. 仅当 $a(u,v)$ 对称正定时两者才一致; 否则, 只能用 Galerkin 法, 而不能用 Ritz 法. Ritz 法的优点是: 力学意义更明显 (尤其是特征值问题), 理论基础比较容易建立.

注5.4当Ritz-Galerkin法用于非齐边值问题时，要根据边值条件的两种不同类型（本质的和自然的）作相应处理.对非齐次自然边值条件，只要适当修改右端即可，不必对基函数加任何限制.对于非齐次本质边值条件，应对它齐次化后再用Ritz-Galerkin方法.例如非齐次边值问题(5.25)，其右端点为自然边值条件，因此右端应改为 $(f,\varphi_j) + p(b)\beta \varphi_j(b)$ （参看5.3.3小节习题2,3).而左端点为本质边值条件，经齐次化后， $u_{n}$ 形如

$$
u _ {n} (x) = u _ {0} (x) + \sum_ {i = 1} ^ {n} c _ {i} \varphi_ {i} (x),
$$

其中 $\varphi_{i}(a) = 0 (i = 1,2,\dots ,n), u_{0}(x)$ 是满足 $u_{0}(a) = \alpha$ 的任一已知函数相应的Ritz-Galerkin方程变成

$$
\sum_ {i = 1} ^ {n} a \left(\varphi_ {i}, \varphi_ {j}\right) c _ {i} = (f, \varphi_ {j}) + p (b) \beta \varphi_ {j} (b) - a \left(u _ {0}, \varphi_ {j}\right), \quad j = 1, 2, \dots , n. \tag {5.66}
$$

实际计算时取 $u_{0}(x) = \alpha \varphi_{0}(x),\varphi_{0}(x)$ 是满足 $\varphi_0(a) = 1$ 的任一函数（为使右端点条件保持不变，要求 $\varphi_0(x)$ 在 $b$ 附近等于0).对于二维边值问题，精确给出 $u_{0}(x)$ 是困难的，一般只能用插值法得到 $u_{0}(x)$ 的近似式.

注5.5 我们曾经指出，Ritz法只能用于解对称正定微分算子方程，而Galerkin法则可解更一般的微分方程。例如两点边值问题（参看定理5.9）：

$$
\left\{ \begin{array}{l} L u = - \frac {\mathrm {d}}{\mathrm {d} x} \left(p \frac {\mathrm {d} u}{\mathrm {d} x}\right) + r \frac {\mathrm {d} u}{\mathrm {d} x} + q u = f, \quad a <   x <   b, \\ u (a) = 0, u ^ {\prime} (b) = 0, \end{array} \right. \tag {5.67}
$$

其中 $p\in C^1 (\bar{I}),p(x)\geqslant p_{\min} > 0,r,q\in C(\bar{I}),f\in L^2 (I),I = (a,b).$ 与之相应的双线性形式为

$$
a (u, v) = \int_ {a} ^ {b} \left(p \frac {\mathrm {d} u}{\mathrm {d} x} \frac {\mathrm {d} v}{\mathrm {d} x} + r \frac {\mathrm {d} u}{\mathrm {d} x} v + q u v\right) \mathrm {d} x.
$$

显然 $a(u,v)$ 非对称正定，除非 $r\equiv 0,q\geqslant 0.$ 因此不能用Ritz法解(5.40).但 Galerkin法仍然可用，且导出的线性方程和（5.65）相同.

下面我们考虑 Galerkin 法的误差估计. 我们有如下证明的 Céa 引理:

引理5.5（Céa）设 $u$ 是变分问题(5.60)的解， $u_{n}$ 是其Galerkin离散(5.62)的

解.假设双线性形式 $a(\cdot ,\cdot)$ 满足强制性(5.57）和连续性(5.58)，则成立如下估计：

$$
\| u - u _ {n} \| _ {1} \leqslant \frac {M}{\gamma} \inf  _ {v \in U _ {n}} \| u - v \| _ {1}. \tag {5.68}
$$

证明 既然 $U_{n} \subset U$ , 由 $u$ 和 $u_{n}$ 的定义得,

$$
a (u, v _ {n}) = (f, v _ {n}), \quad \forall v _ {n} \in U _ {n},
$$

$$
a \left(u _ {n}, v _ {n}\right) = \left(f, v _ {n}\right), \quad \forall v _ {n} \in U _ {n}.
$$

相减得

$$
a \left(u - u _ {n}, v _ {n}\right) = 0, \quad \forall v _ {n} \in U _ {n}, \tag {5.69}
$$

即 Galerkin 解的误差按 $a(\cdot, \cdot)$ 内积与离散空间 $U_{n}$ 垂直，称为“Galerkin 正交性”。从而

$$
\begin{array}{l} \gamma \| u - u _ {n} \| _ {1} ^ {2} \leqslant a (u - u _ {n}, u - u _ {n}) \quad (\text {强 制 性}) \\ = a \left(u - u _ {n}, u - v _ {n}\right) + a \left(u - u _ {n}, v _ {n} - u _ {n}\right) \\ = a \left(u - u _ {n}, u - v _ {n}\right) \quad (\text {正 交 性}) \\ \leqslant M \| u - u _ {n} \| _ {1} \| u - v _ {n} \| _ {1}. \quad (\text {连 续 性}) \\ \end{array}
$$

两边消去共同因子 $\| u - u_n\| _1$ ，并关于 $v_{n}\in U_{n}$ 取下确界即得证明

注5.6 (1) Galekin法的Céa引理并没有要求双线性形式满足对称性，可以应用到非对称问题

(2) Céa 引理指出 Galerkin 方法的解的误差与近似空间 $U_{n}$ 中的最佳逼近的误差同阶. 所以应用 Galerkin 方法的关键是如何构造近似空间 $V_{h}$ 使得其中的函数能够很好的逼近精确解.  
(3) Galerkin 法和相应的 Céa 引理显然可以推广到更一般的边值问题. 事实上, $a(\cdot, \cdot)$ 可以是一般的双线性形式, 不局限于前面那些例子. 另外, 把 $H^{1}$ 空间改为 $H^{m}$ 空间, $H^{1}$ 范数改为 $H^{m}$ 范数, 那么相应修改的 Céa 引理也成立. 比如取 $m = 2$ , 就可以应用到四阶椭圆边值问题, 这里就不展开了.

例5.4 用Ritz-Galerkin法解边值问题

$$
\left\{ \begin{array}{l l} u ^ {\prime \prime} + u = - x, & 0 <   x <   1, \\ u (0) = u (1) = 0. \end{array} \right. \tag {5.70}
$$

此时， $U = H_0^1(I) (I = (0,1))$ ， $H = L^2(I)$ . 于 $H_0^1(I)$ 取一族基函数 $\varphi_i(x) (i = 1,2,\dots)$ ，使每一 $\varphi_i(x)$ 满足齐次边值条件，彼此线性独立，且构成 $H_0^1$ 的完全系统. 以 $\varphi_1, \varphi_2, \dots, \varphi_n$ 为基底张开的子空间就是 $n$ 维空间 $U_n$ .

通常有两种选取 $\varphi_{i}$ 的方法. 一种是选 $\varphi_{i}$ 为三角多项式

$$
\varphi_ {i} (x) = \sin (i \pi x), \quad i = 1, 2, \dots ,
$$

另一种是取 $\varphi_{i}$ 为代数多项式

$$
\varphi_ {i} (x) = \omega (x) x ^ {i - 1}, \quad i = 1, 2, \dots .
$$

为使 $\varphi_{i}$ 满足边值条件，取

$$
\omega (x) = x (1 - x),
$$

将 $u_{n}(x)$ 表示成

$$
u _ {n} (x) = \sum_ {i = 1} ^ {n} c _ {i} \varphi_ {i} (x) = x (1 - x) \left(c _ {1} + c _ {2} x + \dots + c _ {n} x ^ {n - 1}\right).
$$

先令 $n = 1$ ，则 $u_{1} = c_{1}x(1 - x)$ .由（5.65） $(n = 1)$ ， $c_{1}$ 满足方程

$$
c _ {1} \int_ {0} ^ {1} \left(\varphi_ {1} ^ {\prime \prime} + \varphi_ {1}\right) \varphi_ {1} \mathrm {d} x = - \int_ {0} ^ {1} x ^ {2} (1 - x) \mathrm {d} x.
$$

经计算，得

$$
- \frac {3}{1 0} c _ {1} = - \frac {1}{1 2}, \quad c _ {1} = \frac {5}{1 8}, \quad u _ {1} = \frac {5}{1 8} x (1 - x).
$$

再令 $n = 2$ 以 $u_{2} = c_{1}\varphi_{1} + c_{2}\varphi_{2},f = -x$ 代到(5.65)，经简单计算，得Ritz-Galerkin方程：

$$
\left\{ \begin{array}{l} - \frac {3}{1 0} c _ {1} - \frac {3}{2 0} c _ {2} = - \frac {1}{1 2}, \\ - \frac {3}{2 0} c _ {1} - \frac {1 3}{1 0 5} c _ {2} = - \frac {1}{2 0}. \end{array} \right.
$$

解之, 得 $c_{1} = \frac{71}{369}, c_{2} = \frac{7}{41}$ , 于是

$$
u _ {2} = x (1 - x) \left(\frac {7 1}{3 6 9} + \frac {7}{4 1} x\right).
$$

边值问题 (5.70) 的精确解为

$$
u _ {*} = \frac {\sin x}{\sin 1} - x.
$$

表5.1列出 $u_{1}(x), u_{2}(x), u_{*}(x)$ 于 $x = \frac{1}{4}, \frac{1}{2}, \frac{3}{4}$ 的函数值

表5.1 $u_{*}(x),u_{1}(x),u_{2}(x)$ 于 $x = \frac{1}{4},\frac{1}{2},\frac{3}{4}$ 的函数值  

<table><tr><td>x</td><td>u*(x)</td><td>u1(x)</td><td>u2(x)</td></tr><tr><td>1/4</td><td>0.044</td><td>0.052</td><td>0.044</td></tr><tr><td>1/2</td><td>0.070</td><td>0.069</td><td>0.069</td></tr><tr><td>3/4</td><td>0.060</td><td>0.052</td><td>0.060</td></tr></table>

上述例子是简单的. 实际应用中的问题要复杂得多. 例如基函数的选取, 它必须满足本质边值条件. 在有限元方法出现以前, 通常选代数或三角多项式为基函数, 除特别规则的区域外, 要它们满足边值条件是困难的.

下一节我们将对规则区域和周期边值条件介绍一类有效的谱方法和拟谱法. 本书介绍的有限元法, 提供了系统构造基函数或子空间的方法, 可用于求解复杂的边值问题.

# 5.5.1 习题

1. 用Ritz-Galerkin方法求边值问题

$$
\left\{ \begin{array}{l} u ^ {\prime \prime} + u = x ^ {2}, \quad 0 <   x <   1, \\ u (0) = 0, \quad u (1) = 1 \end{array} \right.
$$

的第 $n$ 次近似 $u_{n}(x)$ ，基函数为 $\varphi_i(x) = \sin (i\pi x), i = 1,2,\dots ,n.$

# *5.6 谱方法

本节针对规则区域, 例如一维区间, 二维矩形以及三维长方体等乘积型区域和周期边值条件, 介绍Fourier谱方法, 这是经典Ritz-Galerkin法常用的一种方法. 由于该方法的计算量大, 且要求基函数满足边值条件, 所以在应用中受到很大限制. 1965年, 出现了计算离散Fourier变换的快速算法——FFT算法(参见[6]), 这不仅给Fourier谱方法提供了快速发展的机遇, 而且还将它推广到关于一般正交多项式展开的谱方法(参见[13]). 作为模型, 我们考虑两点边值问题, 推广到高维乘积型区域边值问题并不困难.

# 5.6.1 三角函数逼近

现在假设 $H^{m}(0,2\pi)$ 是定义在 $(0,2\pi)$ 取复值的 Sobolev 空间, $H_{p}^{m}(0,2\pi)$ 是 $H^{m}(0,2\pi)$ 中以 $2\pi$ 为周期的函数组成的 Sobolev 子空间:

$$
H _ {p} ^ {m} (0, 2 \pi) = \left\{f: f \in H ^ {m} (0, 2 \pi), f (x + 2 \pi) = f (x) \right\},
$$

其内积和范数分别为

$$
(f, g) _ {m} = \int_ {0} ^ {2 \pi} \sum_ {s = 0} ^ {m} f ^ {(s)} \bar {g} ^ {(s)} d x
$$

和

$$
\| f \| _ {m} = \sqrt {(f , f) _ {m}}.
$$

往后以 $(f,g)$ 和 $\| f\|$ 分别表示 $(f,g)_0$ 和 $\| f\| _0$

设 $f\in H^{m}(0,2\pi)$ ，将 $f(x)$ 展成Fourier级数：

$$
f (x) = \sum_ {k = - \infty} ^ {\infty} \hat {f} (k) e ^ {i k x}, \tag {5.71}
$$

其中

$$
\hat {f} (k) = \frac {1}{2 \pi} \int_ {0} ^ {2 \pi} f (x) \mathrm {e} ^ {- \mathrm {i} k x} \mathrm {d} x = \frac {1}{2 \pi} \left(f, \mathrm {e} ^ {\mathrm {i} k x}\right)
$$

为 $f$ 的Fourier系数.对(5.71)逐项微商(求广义导数),得

$$
\frac {\mathrm {d} ^ {s} f (x)}{\mathrm {d} x ^ {s}} = \sum_ {k = - \infty} ^ {\infty} \hat {f} (k) (\mathrm {i} k) ^ {s} \mathrm {e} ^ {\mathrm {i} k x}, \quad 0 \leqslant s \leqslant m. \tag {5.72}
$$

由三角函数系 $\left\{\mathrm{e}^{\mathrm{i}kx}\right\}_{k = -\infty}^{\infty}$ 的正交性

$$
\left(\mathrm {e} ^ {\mathrm {i} j x}, \mathrm {e} ^ {\mathrm {i} k x}\right) = \int_ {0} ^ {2 \pi} \mathrm {e} ^ {\mathrm {i} (j - k) x} \mathrm {d} x = \left\{ \begin{array}{l l} 0, & k \neq j, \\ 2 \pi , & k = j, \end{array} \right. \tag {5.73}
$$

可得

$$
\begin{array}{l} \left\| \frac {\mathrm {d} ^ {s} f}{\mathrm {d} x ^ {s}} \right\| ^ {2} = \int_ {0} ^ {2 \pi} \left(\sum_ {j = - \infty} ^ {\infty} \hat {f} (j) (\mathrm {i} j) ^ {s} \mathrm {e} ^ {\mathrm {i} j x}\right) \left(\sum_ {k = - \infty} ^ {\infty} \bar {\hat {f}} (k) (- \mathrm {i} k) ^ {s} \mathrm {e} ^ {- \mathrm {i} k x}\right) \mathrm {d} x \\ = 2 \pi \sum_ {k = - \infty} ^ {\infty} k ^ {2 s} | \hat {f} (k) | ^ {2}. \tag {5.74} \\ \end{array}
$$

取无穷级数 (5.71) 的 $2N + 1$ 项和：

$$
f _ {N} (x) = \sum_ {k = - N} ^ {N} \hat {f} (x) \mathrm {e} ^ {\mathrm {i} k x}, \tag {5.75}
$$

我们自然关心 $f_{N}(x)$ 对 $f(x)$ 的逼近性

引理5.6 设 $f(x) \in H^{m}(0,2\pi)$ ，则对 $s(0 \leqslant s \leqslant m)$ 有估计

$$
\left\| f (x) - f _ {N} (x) \right\| _ {s} \leqslant C N ^ {s - m} \| f \| _ {m},
$$

其中 $C$ 是与 $f, N$ 无关的常数

证明 由(5.71)和(5.75)，

$$
\frac {\mathrm {d} ^ {j} (f (x) - f _ {N} (x))}{\mathrm {d} x ^ {j}} = \frac {\mathrm {d} ^ {j}}{\mathrm {d} x ^ {j}} \sum_ {| n | > N} \hat {f} (n) \mathrm {e} ^ {\mathrm {i} n x} = \sum_ {| n | > N} \hat {f} (n) (\mathrm {i} n) ^ {j} \mathrm {e} ^ {\mathrm {i} n x}.
$$

对 $0 \leqslant j \leqslant s \leqslant m$ , 由 (5.74) 有

$$
\begin{array}{l} \left\| \frac {\mathrm {d} ^ {j} (f (x) - f _ {N} (x))}{\mathrm {d} x ^ {j}} \right\| ^ {2} = 2 \pi \sum_ {n > N} | \hat {f} (n) | ^ {2} n ^ {2 j} = 2 \pi \sum_ {n > N} | \hat {f} (n) | ^ {2} n ^ {2 (j - m)} n ^ {2 m} \\ \leqslant 2 \pi N ^ {2 (j - m)} \sum_ {| n | > N} | \hat {f} (n) | ^ {2} n ^ {2 m} \leqslant N ^ {2 (j - m)} \left\| \frac {\mathrm {d} ^ {m} f}{\mathrm {d} x ^ {m}} \right\| ^ {2}. \\ \end{array}
$$

关于 $j(0 \leqslant j \leqslant s)$ 取和, 两边开方, 即知结论成立.

记 $U_{N} = \operatorname{span}\left\{\mathrm{e}^{\mathrm{i}kx}\right\}_{k = -N}^{N}, P_{N}$ 是由 $H = H^{0}$ 到 $U_{N}$ 的 $L^{2}$ 投影算子, 即 $\forall f \in H$ 有唯一 $P_{N}f \in U_{N}$ , 使

$$
(P _ {N} f, v) = (f, v), \quad \forall v \in U _ {N}.
$$

则由正交性(5.73)知， $f_{N}(x)$ 就是 $f(x)$ 从 $H(0,2\pi)$ 到 $U_{N}$ 的 $L^2$ 投影，记为 $f_{N} = P_{N}f.$ 为应用FFT算法，最好取 $N = 2^{m}:U_{N} = \operatorname {span}\left\{\mathrm{e}^{\mathrm{i}kx}\right\}_{k = -N}^{N - 1}.$ 下面考虑 $f(x)$ 在 $U_{N}$ 的插值逼近.在 $[0,2\pi ]$ 中引入 $2N$ 个等距节点 $x_{m} = \frac{m\pi}{N},m = 0,1,\dots ,2N.$ 利用三角函数的性质可以证明

$$
\frac {1}{2 N} \sum_ {k = - N} ^ {N - 1} \mathrm {e} ^ {\mathrm {i} k x _ {m}} = \left\{ \begin{array}{l l} 1, & m = 0, \\ 0, & m \neq 0. \end{array} \right. \tag {5.76}
$$

令

$$
l _ {m} (x) = \frac {1}{2 N} \sum_ {k = - N} ^ {N - 1} \mathrm {e} ^ {\mathrm {i} k (x - x _ {m})}. \tag {5.77}
$$

由(5.76)有

$$
l _ {m} \left(x _ {n}\right) = \delta_ {m n} = \left\{ \begin{array}{l l} 1, & m = n, \\ 0, & m \neq n. \end{array} \right.
$$

可见 $l_{m}(x)$ 可作为Lagrange插值基函数.这样， $f(x)$ 在 $U_{N}$ 中以 $\{x_j\}$ 为插值节点的插值多项式为

$$
\begin{array}{l} I _ {N} f (x) = \sum_ {m = 1} ^ {2 N} f \left(x _ {m}\right) l _ {m} (x) = \sum_ {m = 1} ^ {2 N} f \left(x _ {m}\right) \frac {1}{2 N} \sum_ {k = - N} ^ {N - 1} \mathrm {e} ^ {\mathrm {i} k \left(x - x _ {m}\right)} \\ = \sum_ {k = - N} ^ {N - 1} \mathrm {e} ^ {\mathrm {i} k x} \frac {1}{2 N} \sum_ {m = 1} ^ {2 N} f \left(x _ {m}\right) \mathrm {e} ^ {- \mathrm {i} k x _ {m}}. \tag {5.78} \\ \end{array}
$$

定义离散内积和范数：

$$
(u, v) _ {N} = \frac {\pi}{N} \sum_ {m = 1} ^ {2 N} u (x _ {m}) \bar {v} (x _ {m}), \| u \| _ {N} ^ {2} = (u, u) _ {N}. \tag {5.79}
$$

则

$$
I _ {N} f (x) = \sum_ {k = - N} ^ {N - 1} \tilde {f} (k) \mathrm {e} ^ {\mathrm {i} k x}, \tag {5.80}
$$

其中

$$
\tilde {f} (k) = \frac {1}{2 \pi} \left(f, \mathrm {e} ^ {\mathrm {i} k x}\right) _ {N} \tag {5.81}
$$

是离散Fourier系数

# 5.6.2 Fourier 谱方法

现在举例说明谱方法的应用

考虑求解周期边值问题：

$$
L u \equiv - u ^ {\prime \prime} + \lambda u = f (x), \quad x \in (0, 2 \pi), \tag {5.82a}
$$

$$
u (0) = u (2 \pi) = 0, \tag {5.82b}
$$

其中 $\lambda > 0$ 是常数, $f(x)$ 为 $2\pi$ 周期函数. 引进空间

$$
U = \left\{u \in H ^ {1} (0, 2 \pi): u (0) = u (2 \pi) = 0 \right\}, \tag {5.83}
$$

用 $v$ 的复共轭函数 $\bar{v} (x)\in U$ 乘方程(5.82a）两端，在 $[0,2\pi ]$ 上积分，并施行分部积分.利用周期性边值条件(5.82b)，可得

$$
\begin{array}{l} \int_ {0} ^ {2 \pi} L u \cdot \bar {v} \mathrm {d} x = \int_ {0} ^ {2 \pi} u ^ {\prime} \bar {v} ^ {\prime} \mathrm {d} x - \left. u ^ {\prime} \bar {v} \right| _ {0} ^ {2 \pi} + \int_ {0} ^ {2 \pi} \lambda u \bar {v} \mathrm {d} x \\ = \int_ {0} ^ {2 \pi} u ^ {\prime} \bar {v} ^ {\prime} d x + \int_ {0} ^ {2 \pi} \lambda u \bar {v} d x = \int_ {0} ^ {2 \pi} f \bar {v} d x. \\ \end{array}
$$

令

$$
a (u, v) = \int_ {0} ^ {2 \pi} \left(u ^ {\prime} \bar {v} ^ {\prime} + \lambda u \bar {v}\right) d x, \quad (f, v) = \int_ {0} ^ {2 \pi} f \bar {v} d x.
$$

于是问题 (5.82a) (5.82b) 的变分形式为: 求 $u \in U$ , 使得

$$
a (u, v) = (f, v), \quad \forall v \in U. \tag {5.84}
$$

在 Galerkin 法中, 取子空间 $U_{n} = \operatorname{span}\left\{\mathrm{e}^{\mathrm{i}kx}\right\}_{k = -N}^{N} \subset U$ , 就导出所谓的 Fourier 谱方法.

现在介绍解边值问题（5.82a）（5.82b）的谱方法．取基函数 $\varphi_{k} = \mathrm{e}^{\mathrm{i}kx} - 1$ ， $k = \pm 1,\pm 2,\dots ,\pm N,U_N = \operatorname {span}\{\varphi_k:k = \pm 1,\pm 2,\dots ,\pm N\}$ ，则谱方法为：求

$$
\begin{array}{l} u _ {N} = \sum_ {\substack {k = - N \\ k \neq 0}} ^ {N} c _ {k} \varphi_ {k} = \sum_ {\substack {k = - N \\ k \neq 0}} ^ {N} c _ {k} \left(\mathrm {e} ^ {\mathrm {i} k x} - 1\right) = \sum_ {k = - N} ^ {N} c _ {k} \mathrm {e} ^ {\mathrm {i} k x}, \tag{5.85} \\ c _ {0} = - \sum_ {k = - N, k \neq 0} ^ {N} c _ {k}, \\ \end{array}
$$

满足

$$
\begin{array}{l} a \left(u _ {N}, \varphi_ {j}\right) = \sum_ {k = - N} ^ {N} c _ {k} a \left(\mathrm {e} ^ {\mathrm {i} k x}, \varphi_ {j}\right) \\ = \sum_ {k = - N} ^ {N} c _ {k} a \left(\mathrm {e} ^ {\mathrm {i} k x}, \mathrm {e} ^ {\mathrm {i} j x}\right) - \sum_ {k = - N} ^ {N} c _ {k} a \left(\mathrm {e} ^ {\mathrm {i} k x}, 1\right), \quad j = \pm 1, \pm 2, \dots , \pm N, \\ = (f, \varphi_ {j}) \tag {5.86} \\ \end{array}
$$

和

$$
\sum_ {k = - N} ^ {N} c _ {k} = 0. \tag {5.87}
$$

由正交性(5.73)得

$$
a \left(\mathrm {e} ^ {\mathrm {i} k x}, \mathrm {e} ^ {\mathrm {i} j x}\right) = \int_ {0} ^ {2 \pi} \left[ - (\mathrm {i} k) (\mathrm {i} j) \mathrm {e} ^ {\mathrm {i} k x} \mathrm {e} ^ {- \mathrm {i} j x} + \lambda \mathrm {e} ^ {\mathrm {i} k x} \mathrm {e} ^ {- \mathrm {i} j x} \right] \mathrm {d} x = \left\{ \begin{array}{l l} 0, & k \neq j, \\ 2 \pi (\lambda + j ^ {2}), & k = j, \end{array} \right. \tag {5.88}
$$

特别地, 取 $k = j = 0$ , 得

$$
a (1, 1) = 2 \pi \lambda .
$$

又

$$
(f, \varphi_ {j}) = \int_ {0} ^ {2 \pi} f (x) \mathrm {e} ^ {- \mathrm {i} j x} \mathrm {d} x - \int_ {0} ^ {2 \pi} f (x) \mathrm {d} x. \tag {5.89}
$$

于是方程组（5.86）化为

$$
2 \pi (\lambda + j ^ {2}) c _ {j} - 2 \pi \lambda c _ {0} = (f, \varphi_ {j}), \quad j = \pm 1, \pm 2, \dots , \pm N.
$$

与(5.87)联立得方程组

$$
\left[ \begin{array}{c c c c c c} 1 & & & & & a _ {- 1} \\ & 1 & & & & a _ {1} \\ & & \ddots & & & \vdots \\ & & & 1 & & a _ {- N} \\ & & & & 1 & a _ {N} \\ - 1 & - 1 & - 1 & - 1 & - 1 \end{array} \right] \left[ \begin{array}{c} c _ {- 1} \\ c _ {1} \\ \vdots \\ c _ {- N} \\ c _ {N} \\ c _ {0} \end{array} \right] = \left[ \begin{array}{c} b _ {- 1} \\ b _ {1} \\ \vdots \\ b _ {- N} \\ b _ {N} \\ 0 \end{array} \right], \tag {5.90}
$$

其中

$$
a _ {j} = - \frac {\lambda}{\lambda + j ^ {2}}, \quad b _ {j} = \frac {(f , \varphi_ {j})}{2 \pi (\lambda + j ^ {2})}, \quad j = \pm 1, \pm 2, \dots , \pm N. \tag {5.91}
$$

用消元法解出系数 $c_{j}$ 即得近似解

由引理5.5和引理5.6得

$$
\| u - u _ {N} \| _ {1} \leqslant C \inf  _ {v _ {n} \in U _ {n}} \| u - v _ {n} \| _ {1} \leqslant C N ^ {- (s - 1)} \| u \| _ {s}.
$$

![](images/2cae9e56879913f0d9792cae372eb052b140c92845724d93a8d208e43b4226ae.jpg)

人物简介

由上式可看出， $u(x)$ 的光滑性越好， $u_{N}(x)$ 收敛得越快.特别地，若 $u(x)$ 是无穷次可微的周期函数，则 $u_{N}(x)$ 收敛于 $u(x)$ 的速度快于 $\frac{1}{N}$ 的任何有限次幂.所以也说谱方法具有“指数收敛性”.

注5.7Fourier谱方法要求解具有周期性.对非周期情形，可作周期性扩张，但在边界点会出现间断，将产生不应有的振荡.为此，人们研究用Chebyshev,Legendre等多项式作为逼近工具，并简称相应的Galerkin方法为谱方法.

例5.5 用谱方法求解

$$
- u ^ {\prime \prime} + u = 2 x \sin x - 2 \cos x, \quad x \in (0, 2 \pi)
$$

$$
u (0) = u (2 \pi) = 0
$$

(精确解 $u = x\sin x$

令 $f(x) = 2x\sin x - 2\cos x$ ，谱方法右端

$$
\begin{array}{l} (f, \varphi_ {j}) = \int_ {0} ^ {2 \pi} f \bar {\varphi} _ {j} \mathrm {d} x = \int_ {0} ^ {2 \pi} f (x) \mathrm {e} ^ {- \mathrm {i} j x} \mathrm {d} x - \int_ {0} ^ {2 \pi} f (x) \mathrm {d} x \\ = \int_ {0} ^ {2 \pi} f (x) \cos j x \mathrm {d} x - \int_ {0} ^ {2 \pi} f (x) \mathrm {d} x - \mathrm {i} \int_ {0} ^ {2 \pi} f (x) \sin j x \mathrm {d} x. \tag {5.92} \\ \end{array}
$$

经计算得

$$
(f, \varphi_ {j}) = \left\{ \begin{array}{l l} 4 \pi \frac {j ^ {2}}{j ^ {2} - 1}, & j \neq \pm 1, \\ \pi - 2 \mathrm {i} \pi^ {2}, & j = 1, \\ \pi + 2 \mathrm {i} \pi^ {2}, & j = - 1. \end{array} \right. \tag {5.93}
$$

另一方面，由(5.91)有（取 $\lambda = 1$ ）

$$
b _ {j} = \frac {(f , \varphi_ {j})}{2 \pi (1 + j ^ {2})} = \left\{ \begin{array}{l l} - \frac {2 j ^ {2}}{(1 - j ^ {4})}, & j \neq \pm 1, \\ \frac {1}{4} - i \frac {\pi}{2}, & j = 1, \\ \frac {1}{4} + i \frac {\pi}{2}, & j = - 1. \end{array} \right. \tag {5.94}
$$

$$
a _ {j} = \frac {- 1}{1 + j ^ {2}}, \quad j = \pm 1, \pm 2, \dots , \pm N. \tag {5.95}
$$

最后解方程组 (5.90). 用消元法得

$$
\left[ \begin{array}{c c c c c} 1 & & & & a _ {- 1} \\ 1 & & & & a _ {1} \\ & \ddots & & & \vdots \\ & & 1 & & a _ {- N} \\ & & & 1 & a _ {N} \\ & & & & d _ {0} \end{array} \right] \left[ \begin{array}{c} c _ {- 1} \\ c _ {1} \\ \vdots \\ c _ {- N} \\ c _ {N} \\ c _ {0} \end{array} \right] = \left[ \begin{array}{c} b _ {- 1} \\ b _ {1} \\ \vdots \\ b _ {- N} \\ b _ {N} \\ b _ {0} \end{array} \right], \tag {5.96}
$$

其中

$$
b _ {0} = \sum_ {k = - N, k \neq 0} ^ {N} b _ {k} = - 4 \sum_ {j = 2} ^ {N} \frac {j ^ {2}}{\left(1 - j ^ {4}\right)} + \frac {1}{2}, \tag {5.97a}
$$

$$
d _ {0} = - 1 + \sum_ {k = - N, k \neq 0} ^ {N} a _ {k} = - \left(1 + 2 \sum_ {k = 1} ^ {N} \frac {1}{1 + k ^ {2}}\right). \tag {5.97b}
$$

由(5.96)回代得

$$
c _ {0} = \frac {b _ {0}}{d _ {0}}, \quad c _ {j} = b _ {j} - a _ {j} c _ {0}. \tag {5.98}
$$

从而

$$
\begin{array}{l} u _ {N} = c _ {0} + \sum_ {j = - N, j \neq 0} ^ {N} c _ {j} \mathrm {e} ^ {\mathrm {i} j x} = c _ {0} + \sum_ {j = - N, j \neq 0} ^ {N} b _ {j} \mathrm {e} ^ {\mathrm {i} j x} - c _ {0} \sum_ {j = - N, j \neq 0} ^ {N} a _ {j} \mathrm {e} ^ {\mathrm {i} j x} \\ = c _ {0} \left(1 - \sum_ {j = - N, j \neq 0} ^ {N} a _ {j} \mathrm {e} ^ {\mathrm {i} j x}\right) + \sum_ {j = - N, j \neq 0} ^ {N} b _ {j} \mathrm {e} ^ {\mathrm {i} j x}, \tag {5.99} \\ \end{array}
$$

其中（利用(5.94)和(5.95))

$$
\sum_ {j = - N, j \neq 0} ^ {N} a _ {j} \mathrm {e} ^ {\mathrm {i} j x} = \sum_ {j = 1} ^ {N} \frac {- 2}{1 + j ^ {2}} \cos j x, \tag {5.100a}
$$

$$
\sum_ {j = - N, j \neq 0} ^ {N} b _ {j} \mathrm {e} ^ {\mathrm {i} j x} = - 4 \sum_ {j = 2} ^ {N} \frac {j ^ {2}}{\left(1 - j ^ {4}\right)} \cos j x + \frac {1}{2} \cos x + \pi \sin x. \tag {5.100b}
$$

将(5.97a)—(5.98)和(5.100a)(5.100b)代入(5.99)，则得近似解：

$$
\begin{array}{l} u _ {N} = c _ {0} \left(1 - \sum_ {j = - N, j \neq 0} ^ {N} a _ {j} \mathrm {e} ^ {\mathrm {i} j x}\right) + \sum_ {j = - N, j \neq 0} ^ {N} b _ {j} \mathrm {e} ^ {\mathrm {i} j x} \\ = \left\{- \left(1 + 2 \sum_ {k = 1} ^ {N} \frac {1}{1 + k ^ {2}}\right) ^ {- 1} \left[ - 4 \sum_ {j = 2} ^ {N} \frac {j ^ {2}}{(1 - j ^ {4})} + \frac {1}{2} \right] \right\}. \\ \left\{1 + 2 \sum_ {j = 1} ^ {N} \frac {1}{1 + j ^ {2}} \cos j x \right\} - 4 \sum_ {j = 2} ^ {N} \frac {j ^ {2}}{\left(1 - j ^ {4}\right)} \cos j x + \frac {1}{2} \cos x + \pi \sin x. \tag {5.101} \\ \end{array}
$$

计算结果如图5.2、图5.3所示，可见谱方法的精度是很高的。

![](images/019c56b986aeb2c240e5ff96762073f3430e569fd0ab591e0f7c4f75bf45dce1.jpg)  
图5.2 $N = 4$ 计算结果

![](images/f9ca9ee11226958d23d10dbb61ec9551bfdeddb1e715ea100ca42a687342054f.jpg)  
图5.3 $N = 8$ 计算结果

# 5.6.3 拟谱方法（配置法）

谱方法要计算许多诸如 $a\left(\mathrm{e}^{\mathrm{i}kx}, \mathrm{e}^{\mathrm{i}jx}\right)$ 的内积, 对变系数方程, 计算量较大, 有时要用数值积分公式. 现在采用配置法, 称为拟谱方法, 可明显减少计算量.

设边值问题为求 $2\pi$ 周期函数 $u$ ，满足

$$
L u = f. \tag {5.102}
$$

拟谱方法是选定节点组 $\{x_{j}\}$ ，求 $u_{N}\in U_{N}$ ，使得 $u_{N}$ 在 $\{x_j\}$ 上满足方程(5.102).

仍以边值问题(5.82a)（5.82b）为例介绍拟谱方法的应用

取基函数 $\varphi_{k} = \mathrm{e}^{\mathrm{i}kx} - 1$ ( $k = \pm 1, \pm 2, \dots, \pm N$ ), $U_{N} = \operatorname{span}\{\varphi_{k}, k = \pm 1, \pm 2, \dots, \pm N\}$ ,

$$
u _ {N} = \sum_ {\substack {k = - N \\ k \neq 0}} ^ {N} c _ {k} \varphi_ {k} = \sum_ {\substack {k = - N \\ k \neq 0}} ^ {N} c _ {k} \left(\mathrm {e} ^ {\mathrm {i} k x} - 1\right) = \sum_ {k = - N} ^ {N} c _ {k} \mathrm {e} ^ {\mathrm {i} k x}, \tag{5.103}
$$

$$
c _ {0} = - \sum_ {k = - N, k \neq 0} ^ {N} c _ {k}.
$$

则

$$
\frac {\mathrm {d} u _ {N}}{\mathrm {d} x} = \sum_ {k = - N} ^ {N} c _ {k} (\mathrm {i} k) \mathrm {e} ^ {\mathrm {i} k x}, \quad \frac {\mathrm {d} ^ {2} u _ {N}}{\mathrm {d} x ^ {2}} = \sum_ {k = - N} ^ {N} c _ {k} (- k ^ {2}) \mathrm {e} ^ {\mathrm {i} k x}.
$$

由于 $\mathrm{e}^{\mathrm{i}kx_0} = \mathrm{e}^{\mathrm{i}kx_{2N}}$ ，故可设配置点为 $x_{m} = \frac{m\pi}{N}, m = 1,2,\dots ,2N,$ 令 $u_{N}$ 在 $2N$ 个配置点 $x_{m}$ 上满足方程（5.82a）(5.82b)，得到

$$
\sum_ {k = - N} ^ {N} c _ {k} k ^ {2} \mathrm {e} ^ {\mathrm {i} k x _ {m}} + \lambda \sum_ {k = - N} ^ {N} c _ {k} \mathrm {e} ^ {\mathrm {i} k x _ {m}} = f (x _ {m}), \quad m = 1, 2, \dots , 2 N. \tag {5.104}
$$

两端乘 $\mathrm{e}^{-\mathrm{i}jx_m}$ , 并关于 $m = 1,2,\dots ,2N$ 求和, 则左端为

$$
\begin{array}{l} \sum_ {m = 1} ^ {2 N} \sum_ {k = - N} ^ {N} c _ {k} k ^ {2} \mathrm {e} ^ {\mathrm {i} k x _ {m}} \mathrm {e} ^ {- \mathrm {i} j x _ {m}} + \lambda \sum_ {m = 1} ^ {2 N} \sum_ {k = - N} ^ {N} c _ {k} \mathrm {e} ^ {\mathrm {i} k x _ {m}} \mathrm {e} ^ {- \mathrm {i} j x _ {m}} \\ = \sum_ {k = - N} ^ {N} c _ {k} k ^ {2} \sum_ {m = 1} ^ {2 N} \mathrm {e} ^ {\mathrm {i} k x _ {m}} \mathrm {e} ^ {- \mathrm {i} j x _ {m}} + \lambda \sum_ {k = - N} ^ {N} c _ {k} \sum_ {m = 1} ^ {2 N} \mathrm {e} ^ {\mathrm {i} k x _ {m}} \mathrm {e} ^ {- \mathrm {i} j x _ {m}} \\ = \sum_ {k = - N} ^ {N} c _ {k} k ^ {2} \frac {N}{\pi} \left(\mathrm {e} ^ {\mathrm {i} k x}, \mathrm {e} ^ {\mathrm {i} j x}\right) _ {N} + \lambda \sum_ {k = - N} ^ {N} c _ {k} \frac {N}{\pi} \left(\mathrm {e} ^ {\mathrm {i} k x}, \mathrm {e} ^ {\mathrm {i} j x}\right) _ {N}, \\ \end{array}
$$

右端为

$$
\sum_ {m = 1} ^ {2 N} f \left(x _ {m}\right) \mathrm {e} ^ {- \mathrm {i} j x _ {m}} = \frac {N}{\pi} \left(f, \mathrm {e} ^ {\mathrm {i} j x}\right) _ {N},
$$

其中

$$
(u, v) _ {N} = \frac {\pi}{N} \sum_ {m = 1} ^ {2 N} u (x _ {m}) \bar {v} (x _ {m}).
$$

于是方程组 (5.104) 化为

$$
\sum_ {\substack {k = - N \\ k \neq 0}} ^ {N} c _ {k} \left(k ^ {2} + \lambda\right) \left(\mathrm {e} ^ {\mathrm {i} k x}, \mathrm {e} ^ {\mathrm {i} j x}\right) _ {N} = \left(f, \mathrm {e} ^ {\mathrm {i} j x}\right) _ {N}, \quad j = - N, - N + 1, \dots , N \tag{5.105}
$$

与

$$
c _ {0} = - \sum_ {k = - N, k \neq 0} ^ {N} c _ {k}.
$$

联立求出 $c_{k}(k = 0,\pm 1,\dots ,\pm N)$ ，代到(5.103）即得 $u_{N}$

(5.105) 是 Fourier 谱方法中方程 (5.86) 的离散形式. 计算形如 $\left(\mathrm{e}^{\mathrm{i}kx}, \mathrm{e}^{\mathrm{i}jx}\right)_N$ , $\left(f, \mathrm{e}^{\mathrm{i}jx}\right)_N$ 的离散内积时, 可采用 FFT 算法, 参见 [6]. FFT 算法可将要计算的复数运算由 $O(N^2)$ 个减少到 $O(N \cdot \log_2 N)$ 个. 可以证明, 拟谱方法和谱方法有同样的收敛阶. 对于具复杂系数的方程, 特别是非线性问题, 拟谱方法更为实用 [13].

![](images/ca020d26efbc4ab33eb8a6a2e8e50af53baaaaf3dc2c8c6ab9cfb29903bfa587.jpg)

# 第6章

# 有限元法

有限元法, 实质上就是一种特殊的 Ritz-Galerkin 法, 它和传统的 Ritz-Galerkin 法的主要区别在于, 它用样条函数方法提供了一种选取“局部基函数”或“分片多项式空间”的技术, 从而在很大程度上克服了 Ritz-Galerkin 法选取基函数的困难. 有限元法首先成功地用于结构力学和固体力学, 后又用于流体力学、物理学和其他工程科学 [10, 18, 24, 29, 31, 35]. 现在, 有限元法和差分法一样, 已成为求解偏微分方程, 特别是椭圆型偏微分方程的一种有效数值方法.

有限元法的基本问题可归纳为

(1) 把问题转化成变分形式  
(2) 选定单元的形状, 对求解域作剖分.

一维情形的单元是小区间. 二维情形的重要单元有两种: 四边形 (矩形、任意凸四边形) 和三角形. 三维单元就更复杂多样了, 比如四面体元, 六面体元等. 本书只讨论一维、二维和三维单元.

(3) 构造基函数或单元形状函数  
(4) 形成有限元方程 (Ritz-Galerkin 方程).  
(5) 提供有限元方程的有效解法  
(6) 收敛性及误差估计.

第一个问题已在第5章讲过了，问题(5)见数值线性代数的相关内容，可以采用学过的直接法或迭代法求解.也可以考虑第10章将要介绍的多重网格计算技术.为了便于读者理解，我们先讲一维域上两点边值问题的有限元法，然后推广到二维和三维区域上的二阶椭圆型边值问题.最后，简单介绍如何将有限元法推广到初边值问题，包括抛物型方程和二阶双曲型方程

另外, 为了记号简单和二、三维叙述的统一, 本章黑体的坐标向量 $\pmb{x} \in \mathbb{R}^d$ , 在二维时表示 $\pmb{x} = (x, y)$ , 三维时表示 $\pmb{x} = (x, y, z)$ . 简记 $\| \cdot \|_1 = \| \cdot \|_{H^1(\Omega)}, \| \cdot \| = \| \cdot \|_{L^2(\Omega)}$ .

# 6.1 一维例子

# 6.1.1 两点边值问题及其变分公式

考虑两点边值问题：给定函数 $f(x)$ ，求 $u(x)$ 使得

$$
\left\{ \begin{array}{l} - u ^ {\prime \prime} + u = f (x), \quad 0 <   x <   1, \\ u (0) = 0, \quad u ^ {\prime} (1) = 0. \end{array} \right. \tag {6.1}
$$

如果 $u$ 是(6.1)的解， $v(x)$ 充分光滑且满足 $v(0) = 0$ ，那么由分部积分得

$$
\begin{array}{l} \int_ {0} ^ {1} (- u ^ {\prime \prime} v + u v) d x = - u ^ {\prime} (1) v (1) + u ^ {\prime} (0) v (0) + \int_ {0} ^ {1} \left(u ^ {\prime} v ^ {\prime} + u v\right) d x \\ = \int_ {0} ^ {1} f v \mathrm {d} x. \\ \end{array}
$$

定义双线性形式

$$
a (u, v) = \int_ {0} ^ {1} \left(u ^ {\prime} v ^ {\prime} + u v\right) \mathrm {d} x,
$$

及空间

$$
V = H _ {E} ^ {1} = \left\{v \in H ^ {1} (\Omega): v (0) = 0 \right\},
$$

其中 $\Omega = (0,1)$ ，则得(6.1）的变分公式：求 $u\in V$ 使得

$$
a (u, v) = (f, v), \quad \forall v \in V, \tag {6.2}
$$

其中

$$
(f, v) = \int_ {0} ^ {1} f (x) v (x) d x.
$$

由定理5.8知：若 $f\in C([0,1])$ 且 $u\in C^2 ([0,1])$ ，则原问题(6.1）与变分问题(6.2）等价

# 6.1.2 有限元方法

首先将区间 $[0,1]$ 分成 $n$ 份：

$$
0 = x _ {0} <   x _ {1} <   x _ {2} <   \dots <   x _ {n - 1} <   x _ {n} = 1.
$$

记 $K_{i} = [x_{i - 1},x_{i}]$ 为第 $_i$ 个小区间， $h_i = x_i - x_{i - 1}$ 为其长度.定义 $h = \max_{1\leqslant i\leqslant n}h_i$ 称 $\mathcal{M}_h = \{K_i,i = 1,2,\dots ,n\}$ 为 $\varOmega$ 的一个剖分.

我们用 $\mathcal{M}_h$ 上的连续的分段线性函数逼近 $u(x)$ 。引入有限元空间

$$
U _ {h} = \left\{v _ {h} \in C (\Omega): v _ {h} | _ {K _ {i}} \text {是 线 性 多 项 式 ，} i = 1, 2, \dots , n, v _ {h} (0) = 0 \right\}. \tag {6.3}
$$

显然 $U_{h}\subset V$

相应于 (6.2) 的有限元方法为: 求 $u_{h} \in U_{h}$ , 使得

$$
a \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in U _ {h}. \tag {6.4}
$$

可以看出，有限元方法就是 Galerkin 方法的一种（见第 5.5 节），只不过试探函数空间取成了分段线性函数空间。

# 6.1.3 有限元方程组

先引入空间 $U_{h}$ 的一组基函数，称为节点基：对 $i = 1,2,\dots ,n$ ，定义 $\phi_i\in U_h$ 满足

$$
\phi_ {i} \left(x _ {j}\right) = \delta_ {i j} = \left\{ \begin{array}{l l} 1, & i = j, \\ 0, & i \neq j, \end{array} \right.
$$

如图6.1，有如下表达式：

$$
\phi_ {i} = \left\{ \begin{array}{l l} { \frac {x - x _ {i - 1}}{h _ {i}},} & {x _ {i - 1} \leqslant x \leqslant x _ {i},} \\ { \frac {x _ {i + 1} - x}{h _ {i + 1}},} & {x _ {i} <   x \leqslant x _ {i + 1}, \qquad 1 \leqslant i \leqslant n - 1.} \\ {0,} & {x <   x _ {i - 1} \text {或} x > x _ {i + 1},} \end{array} \right.
$$

$$
\phi_ {n} = \left\{ \begin{array}{l l} \frac {x - x _ {n - 1}}{h _ {n}}, & x _ {n - 1} \leqslant x \leqslant 1, \\ 0, & x <   x _ {n - 1}. \end{array} \right.
$$

对任一 $v_{h}\in U_{h}$ ，记 $\boldsymbol{v}_{i}$ 为 $v_{h}$ 在节点 $x_{i}$ 的值，即 $v_{i} = v_{h}(x_{i}),i = 1,2,\dots ,n,$ 则显然

$$
v _ {h} = v _ {1} \phi_ {1} (x) + v _ {2} \phi_ {2} (x) + \dots + v _ {n} \phi_ {n} (x).
$$

![](images/da76123ce56f38751d144d6baeb5d6e6749a3e550dacc31a0e22ff3238232f40.jpg)

![](images/0bf05fcbcf03fb7e05ba337713ab445227a2ebf745c049c2c948d94809d9a7de.jpg)  
图6.1 节点基函数 $\phi_{i}$

我们有

$$
u _ {h} = u _ {1} \phi_ {1} + u _ {2} \phi_ {2} + \dots + u _ {n} \phi_ {n}, \quad \text {其 中} \quad u _ {i} = u _ {h} \left(x _ {i}\right).
$$

在 (6.4) 中取 $v_{h} = \phi_{i}(i = 1,2,\dots ,n)$ ，得关于未知数 $u_{1},u_{2},\dots ,u_{n}$ 的方程组：

$$
a \left(\phi_ {1}, \phi_ {i}\right) u _ {1} + a \left(\phi_ {2}, \phi_ {i}\right) u _ {2} + \dots + a \left(\phi_ {n}, \phi_ {i}\right) u _ {n} = \left(f, \phi_ {i}\right), \quad i = 1, 2, \dots , n. \tag {6.5}
$$

记

$$
a _ {i j} = a \left(\phi_ {j}, \phi_ {i}\right) = \int_ {0} ^ {1} \left(\phi_ {j} ^ {\prime} \phi_ {i} ^ {\prime} + \phi_ {j} \phi_ {i}\right) d x, \quad f _ {i} = (f, \phi_ {i}),
$$

及

$$
\boldsymbol {A} = \left(a _ {i j}\right) _ {n \times n}, \quad \boldsymbol {F} = \left(f _ {i}\right) _ {n \times 1}, \quad \boldsymbol {U} = \left(u _ {i}\right) _ {n \times 1},
$$

则(6.5)改写为

$$
A U = F.
$$

这里 $\mathbf{A}$ 称为刚度矩阵. 显然, 当 $x_{i}$ 和 $x_{j}$ 不相邻时, $a(\phi_j,\phi_i) = 0$ . 所以 $\mathbf{A}$ 是稀疏矩阵.

下面考虑 $a(\phi_j,\phi_i)$ 的计算，事实上

$$
a \left(\phi_ {j}, \phi_ {i}\right) = \sum_ {k = 1} ^ {n} \int_ {x _ {k - 1}} ^ {x _ {k}} \left(\phi_ {j} ^ {\prime} \phi_ {i} ^ {\prime} + \phi_ {j} \phi_ {i}\right) d x.
$$

先计算 $[x_{i - 1},x_i]$ 上的单元刚度矩阵 $A^i = \left(a_{lm}^i\right)_{2\times 2}$

$$
\begin{array}{l} a _ {1 1} ^ {i} := \int_ {x _ {i - 1}} ^ {x _ {i}} \left(\phi_ {i - 1} ^ {\prime} \phi_ {i - 1} ^ {\prime} + \phi_ {i - 1} \phi_ {i - 1}\right) d x = \frac {1}{h _ {i}} + \frac {h _ {i}}{3}, \\ a _ {2 2} ^ {i} := \int_ {x _ {i - 1}} ^ {x _ {i}} \left(\phi_ {i} ^ {\prime} \phi_ {i} ^ {\prime} + \phi_ {i} \phi_ {i}\right) d x = \frac {1}{h _ {i}} + \frac {h _ {i}}{3}, \\ a _ {1 2} ^ {i} = a _ {2 1} ^ {i} := \int_ {x _ {i - 1}} ^ {x _ {i}} \left(\phi_ {i} ^ {\prime} \phi_ {i - 1} ^ {\prime} + \phi_ {i} \phi_ {i - 1}\right) d x = - \frac {1}{h _ {i}} + \frac {h _ {i}}{6}. \\ \end{array}
$$

因此

$$
\begin{array}{l} a \left(\phi_ {i}, \phi_ {i}\right) = \left\{ \begin{array}{l} a _ {2 2} ^ {i} + a _ {1 1} ^ {i + 1} = \frac {1}{h _ {i}} + \frac {1}{h _ {i + 1}} + \frac {h _ {i} + h _ {i + 1}}{3}, \quad i = 1, 2, \dots , n - 1, \\ a _ {2 2} ^ {n} = \frac {1}{h _ {n}} + \frac {h _ {n}}{3}, \quad i = n, \end{array} \right. \\ a \left(\phi_ {i}, \phi_ {i - 1}\right) = a _ {1 2} ^ {i} = - \frac {1}{h _ {i}} + \frac {h _ {i}}{6}, \quad i = 2, 3, \dots , n. \\ \end{array}
$$

再由 (6.5) 得有限元方程组

$$
\left\{ \begin{array}{l} \left(\frac {h _ {1} + h _ {2}}{3} + \frac {1}{h _ {1}} + \frac {1}{h _ {2}}\right) u _ {1} + \left(\frac {h _ {2}}{6} - \frac {1}{h _ {2}}\right) u _ {2} = f _ {1}, \\ \left(\frac {h _ {i}}{6} - \frac {1}{h _ {i}}\right) u _ {i - 1} + \left(\frac {h _ {i} + h _ {i + 1}}{3} + \frac {1}{h _ {i}} + \frac {1}{h _ {i + 1}}\right) u _ {i} + \left(\frac {h _ {i + 1}}{6} - \frac {1}{h _ {i + 1}}\right) u _ {i + 1} = f _ {i}, \\ i = 2, 3, \dots , n - 1, \\ \left(\frac {h _ {n}}{6} - \frac {1}{h _ {n}}\right) u _ {n - 1} + \left(\frac {h _ {n}}{3} + \frac {1}{h _ {n}}\right) u _ {n} = f _ {n}. \end{array} \right.
$$

这是一个三对角方程组，可以用追赶法快速求解

# 6.1.4 先验误差估计

先考虑插值误差估计. 给定 $u \in C([0,1])$ ，其有限元插值 $u_{I} \in U_{h}$ 定义为

$$
u _ {I} = \sum_ {i = 1} ^ {n} u (x _ {i}) \phi_ {i}.
$$

显然 $u_{I}(x_{i}) = u(x_{i}), i = 0,1,\dots ,n,$ 且

$$
u _ {I} (x) = \frac {x _ {i} - x}{h _ {i}} u \left(x _ {i - 1}\right) + \frac {x - x _ {i - 1}}{h _ {i}} u \left(x _ {i}\right), \quad x \in \left[ x _ {i - 1}, x _ {i} \right].
$$

定理6.1 存在常数 $C$ 与 $h_i$ 和 $u$ 无关使得下列估计成立：

$$
\left\| u ^ {\prime} - u _ {I} ^ {\prime} \right\| _ {L ^ {2} \left(K _ {i}\right)} \leqslant C h _ {i} \left\| u ^ {\prime \prime} \right\| _ {L ^ {2} \left(K _ {i}\right)}, \tag {6.6}
$$

$$
\left\| u - u _ {I} \right\| _ {L ^ {2} \left(K _ {i}\right)} \leqslant C h _ {i} ^ {2} \left\| u ^ {\prime \prime} \right\| _ {L ^ {2} \left(K _ {i}\right)}. \tag {6.7}
$$

证明 仅证 (6.6). 首先将 (6.6) 化为参考单元 $\hat{K} = [0,1]$ 上的不等式. 从 $\hat{K}$ 到 $K_{i}$ 的变换为 $x = x_{i - 1} + h_{i}\hat{x}$ , 记 $e(x) = u(x) - u_{I}(x), \hat{e} (\hat{x}) = e(x)$ . 不等式 (6.6) 等价于

$$
\left\| \frac {\mathrm {d} \hat {e}}{\mathrm {d} \hat {x}} \right\| _ {L ^ {2} (\hat {K})} ^ {2} = h _ {i} \left\| \frac {\mathrm {d} e}{\mathrm {d} x} \right\| _ {L ^ {2} (K _ {i})} ^ {2} \leqslant C h _ {i} ^ {3} \| u ^ {\prime \prime} \| _ {L ^ {2} (K _ {i})} ^ {2} = C h _ {i} ^ {3} \| e ^ {\prime \prime} \| _ {L ^ {2} (K _ {i})} ^ {2} = C \left\| \frac {\mathrm {d} ^ {2} \hat {e}}{\mathrm {d} \hat {x} ^ {2}} \right\| _ {L ^ {2} (\hat {K})} ^ {2}.
$$

注意到 $\int_0^1\frac{\mathrm{d}\hat{e}}{\mathrm{d}\hat{x}}\mathrm{d}x = 0,$ 由Poincaré不等式知上式成立，证毕

简记 $\| \cdot \| _1 = \| \cdot \|_{H^1 (\varOmega)},\| \cdot \| = \| \cdot \|_{L^2 (\varOmega)}$ .我们有如下整体插值误差估计，其证明略去.

推论6.1 存在常数 $C$ 与 $h_i$ 和 $u$ 无关使得下列估计成立：

$$
\begin{array}{l} \left\| u - u _ {I} \right\| _ {1} \leqslant C h \left\| u ^ {\prime \prime} \right\|, \\ \left\| u - u _ {I} \right\| \leqslant C h ^ {2} \left\| u ^ {\prime \prime} \right\|. \\ \end{array}
$$

显然

$$
a (u, v) \leqslant \| u \| _ {1} \| v \| _ {1}, \quad a (v, v) = \| v \| _ {1} ^ {2}.
$$

故由引理5.5得

$$
\left\| u - u _ {h} \right\| _ {1} \leqslant \inf  _ {v _ {h} \in U _ {h}} \left\| u - v _ {h} \right\| _ {1} \leqslant \left\| u - u _ {I} \right\| _ {1}.
$$

再由推论6.1得如下 $H^1$ 误差估计：

# 定理6.2

$$
\left\| u - u _ {h} \right\| _ {1} \leqslant C h \left\| u ^ {\prime \prime} \right\|.
$$

因为上面的估计依赖于精确解 $u$ ，所以称其为先验误差估计.下面考虑 $L^2$ 范数下的误差估计.引入对偶问题

$$
- w ^ {\prime \prime} + w = u - u _ {h}, \quad 0 <   x <   1, \quad w (0) = 0, \quad w ^ {\prime} (1) = 0.
$$

两边乘 $-w''$ 并在 $\Omega$ 上积分得：

$$
\left\| w ^ {\prime \prime} \right\| ^ {2} + \left\| w ^ {\prime} \right\| ^ {2} = - (u - u _ {h}, w ^ {\prime \prime}) \leqslant \left\| u - u _ {h} \right\| \left\| w ^ {\prime \prime} \right\|.
$$

从而

$$
\left\| w ^ {\prime \prime} \right\| \leqslant \left\| u - u _ {h} \right\|.
$$

我们有

$$
\| u - u _ {h} \| ^ {2} = \int_ {0} ^ {1} (u - u _ {h}) (- w ^ {\prime \prime} + w) d x
$$

$$
\begin{array}{l} = a (u - u _ {h}, w) = a (u - u _ {h}, w - w _ {I}) \\ \leqslant \| u - u _ {h} \| _ {1} \| w - w _ {I} \| _ {1} \\ \leqslant \| u - u _ {h} \| _ {1} C h \| w ^ {\prime \prime} \| \\ \leqslant C h \| u - u _ {h} \| _ {1} \| u - u _ {h} \|, \\ \end{array}
$$

因此

$$
\left\| u - u _ {h} \right\| \leqslant C h \left\| u - u _ {h} \right\| _ {1}.
$$

由定理6.2即得如下 $L^2$ 误差估计

# 定理6.3

$$
\| u - u _ {h} \| \leqslant C h ^ {2} \| u ^ {\prime \prime} \|.
$$

可以看出有限元解的 $H^{1}$ 和 $L^2$ 误差与相应的插值误差同阶.当精确解未知时，其插值也未知，而有限元解是可计算的.当然，两点边值问题(6.1）的精确解是可以写出来的.然而有限元方法及其理论可以推广应用到更一般的问题，其精确解往往是不知道的.接下来的几节，我们将把本节介绍的一维问题的线性有限元方法推广到更一般的情形

# 6.2 有限元空间的构造

有限元空间按如下步骤构造：

- 将区域 $\Omega$ 剖分成有限个小区域 (单元);  
- 在 (每个) 单元上定义有限元函数 (一般取多项式) 和自由度;  
- 把每个单元上定义的有限元函数拼起来形成有限元空间.

下面先介绍如何定义单元上的有限元函数和自由度

# 6.2.1 有限元及有限元空间

我们先以上一节分段线性有限元函数为例, 看看为了构造整体的有限元空间, 如何选取和确定一个单元上的逼近函数, 再抽象出来一般“有限元”的定义.

简单来说, 上一节的每个单元上的有限元函数构造, 就是给定小线段 $K = A_{1}A_{2}$ (比如 $A_{1}, A_{2}$ 分别为第 $i$ 个小区间 $[x_{i-1}, x_{i}]$ 的两个端点), 取其上的逼近函数空间为线性多项式空间 $\mathcal{P} = P_{1}(K)$ . 当然, 任给 $v \in \mathcal{P}$ , 我们还要考虑如何确定它, 由两点确定一条直线, 我们需要两个自由度, 在上一节, 我们取成了 $v$ 在端点处的值, 即 $v(A_{1})$ 和 $v(A_{2})$ . 让我们稍微深入地理解一下这两个自由度: 任给 $v \in \mathcal{P}$ , 得到 $v(A_{i}) \in \mathbb{R}$ , 实际上建立了一

个从 $\mathcal{P}$ 到 $\mathbb{R}$ 的映射, 记为 $N_{i}(v) = v(A_{i}), i = 1,2$ . 显然, $N_{i}$ 是线性的, 即 $N_{i} \in \mathcal{P}'$ ( $\mathcal{P}$ 的对偶空间), 且任意一个 $\mathcal{P}$ 到 $\mathbb{R}$ 的线性映射都可以表示为 $N_{1}$ 和 $N_{2}$ 的线性组合. 也就是说, 为了描述一个单元上的有限元逼近函数, 我们需要三个要素: 单元 $K$ ; 其上的逼近函数空间 $\mathcal{P}$ ; 为了确定 $\mathcal{P}$ 中函数所需的自由度集合.

定义6.1（有限元）有限元是一个三元组 $(K,\mathcal{P},\mathcal{N})$ 满足下列条件：

(i) $K\subset \mathbb{R}^d$ 为具有分片光滑边界的闭区域(称为单元);   
(ii) $\mathcal{P}$ 是 $K$ 上的有限维函数空间, 记 $n = \dim \mathcal{P}$ ( $\mathcal{P}$ 中的元素称为形状函数);  
(iii) $\mathcal{N} = \{N_1, N_2, \dots, N_n\}$ 是对偶空间 $\mathcal{P}'$ 的一组基底 (称为节点变量或自由度).

有限元的定义中 (i) 和 (ii) 都是容易理解和验证的, 下面给出 (iii) 的等价条件.

引理6.1 设 $\mathcal{P}$ 是 $K$ 上的 $n$ 维函数空间, 则 (iii) 与 (iii)' 等价:

(iii)' 若 $v \in \mathcal{P}$ 满足 $N_{i}(v) = 0, i = 1,2,\dots ,n,$ 则 $v \equiv 0$

证明 设 $\{\varphi_1, \varphi_2, \dots, \varphi_n\}$ 是 $\mathcal{P}$ 的一组基，则任一 $v \in \mathcal{P}$ 可以表示为基函数的线性组合：

$$
v = y _ {1} \varphi_ {1} + y _ {2} \varphi_ {2} + \dots + y _ {n} \varphi_ {n}.
$$

显然 $(\mathrm{iii})'$ 等价于下面方程组只有零解：

$$
y _ {1} N _ {i} \left(\varphi_ {1}\right) + y _ {2} N _ {i} \left(\varphi_ {2}\right) + \dots + y _ {n} N _ {i} \left(\varphi_ {n}\right) = 0, \quad i = 1, 2, \dots , n.
$$

等价于系数矩阵 $M = (N_{i}(\varphi_{j}))_{n\times n}$ 可逆

(iii) 等价于: 对任一 $L \in \mathcal{P}^{\prime}$ , 存在数 $z_{1}, z_{2}, \dots, z_{n}$ , 使得 $L = z_{1}N_{1} + z_{2}N_{2} + \dots + z_{n}N_{n}$ , 等价于下面方程组

$$
z _ {1} N _ {1} \left(\varphi_ {i}\right) + z _ {2} N _ {2} \left(\varphi_ {i}\right) + \dots + z _ {n} N _ {n} \left(\varphi_ {i}\right) = L \left(\varphi_ {i}\right), \quad i = 1, 2, \dots , n
$$

有解, 亦等价于其系数矩阵 $M^{\mathrm{T}}$ 可逆. 证毕

有限元三元组中的 $\mathcal{N}$ 是对偶空间 $\mathcal{P}'$ 的一组基. 为了表示有限元函数, 我们还需要 $\mathcal{P}$ 的一组基. 回忆前一节中引入的节点基函数, 限制在线段 $A_{1}A_{2}$ 上只有两个非零, 记为 $\phi_{1}(x), \phi_{2}(x) \in \mathcal{P}$ , 满足: $\phi_{1}(A_{1}) = 1$ , $\phi_{1}(A_{2}) = 0$ ; $\phi_{2}(A_{1}) = 0$ , $\phi_{2}(A_{2}) = 1$ , 用自由度表示就是 $N_{i}(\phi_{j}) = \delta_{ij}, i,j = 1,2$ . 推广到一般有限元, 我们有如下节点基的定义:

定义6.2（节点基）设 $(K,\mathcal{P},\mathcal{N})$ 是一个有限元．设 $\{\phi_1,\phi_2,\dots ,\phi_n\}$ 为 $\mathcal{P}$ 的一组基，且与 $\mathcal{N}$ 对偶，即 $N_{i}(\phi_{j}) = \delta_{ij}$ ，则称其为 $\mathcal{P}$ 的节点基

显然对任意 $v\in \mathcal{P}$ ，有

$$
v (x) = \sum_ {i = 1} ^ {n} N _ {i} (v) \phi_ {i} (x).
$$

下面考虑有限元插值的定义. 回忆前一节的分段线性有限元空间的插值, 任给 $K = A_{1}A_{2}$ 上的连续函数 $v(x)$ , 其插值就是把两个端点 $(A_{1}, v(A_{1})), (A_{2}, v(A_{2}))$ 连起来的线段所对应的线性多项式, 记为 $I_{K}v$ . 显然有

$$
I _ {K} v = v (A _ {1}) \phi_ {1} + v (A _ {2}) \phi_ {2} = N _ {1} (v) \phi_ {1} + N _ {2} (v) \phi_ {2}.
$$

一般地，我们有如下定义：

定义6.3（局部有限元插值）给定有限元 $(K,\mathcal{P},\mathcal{N})$ ，设 $\{\phi_i,1\leqslant i\leqslant n\} \subset \mathcal{P}$ 为 $\mathcal{P}$ 的节点基.如果函数 $v$ 使得 $N_{i}(v)$ 有定义，这里 $N_{i}\in \mathcal{N},i = 1,2,\dots ,n,$ 那么定义其局部有限元插值为

$$
I _ {K} v := \sum_ {i = 1} ^ {n} N _ {i} (v) \phi_ {i}.
$$

显然 $I_{K}$ 是线性算子. $\mathcal{V}$ 与 $I_{K}v$ 的各个自由度相等，即 $N_{i}(v) = N_{i}(I_{K}v),i = 1,2,\dots ,n,$ 且当 $v\in \mathcal{P}$ 时，有 $I_K v = v$

像前一节一维的例子那样，我们将每个单元上的有限元函数拼起来形成有限元空间。为此，先把一维区间网格剖分的概念推广到一般情形：

定义6.4（网格剖分）设 $\Omega$ 是 $\mathbb{R}^d$ 中的有界区域，称 $\mathcal{M}_h = \{K_j\}_{j=1}^J$ 是 $\Omega$ 的一个剖分，如果每个 $K_i$ 都满足有限元定义中的条件 (i); $\bigcup_{j=1}^{J} K_j = \bar{\Omega}$ ；且 $K_i^\circ \cap K_j^\circ = \phi, \forall i \neq j,$ 这里 $K_i^\circ$ 表示单元 $K_i$ 的内部。

对二维情形，如果每个单元 $K_{i}$ 都是三角形，那么称 $\mathcal{M}_h$ 为 $\Omega$ 的一个三角剖分.如果每个单元 $K_{i}$ 都是四边形，那么称 $\mathcal{M}_h$ 为 $\Omega$ 的一个四边形剖分.对三维情形，如果每个单元 $K_{i}$ 都是四面体，那么称 $\mathcal{M}_h$ 为 $\Omega$ 的一个四面体剖分，或与二维情形统称为三角剖分.

给定 $\Omega$ 的一个剖分, 在每个单元上定义有限元, 然后再拼起来, 就可以构造出 $\Omega$ 上的有限元试探函数空间 $U_{h}$ . 下面定理给出了 $U_{h} \subset H^{m}(\Omega)$ 的一个充要条件.

定理6.4 设 $\Omega$ 是 $\mathbb{R}^d$ 中的有界区域, $\mathcal{M}_h = \{K_j\}_{j=1}^J$ 是 $\Omega$ 的一个剖分. 设 $m \geqslant 1$ , 函数 $v$ 满足 $v|_{K_j} \in C^m(K_j)$ , $1 \leqslant j \leqslant J$ , 则 $v \in H^m(\Omega)$ 的充要条件是 $v \in C^{m-1}(\Omega)$ .

证明 仅证 $m = 1$ 的情形. 对 $m > 1$ , 可以对 $m - 1$ 阶导数应用 $m = 1$ 的结论.

先证充分性. 设 $v \in C(\Omega)$ . 定义函数 $\pmb{w}$ 满足

$$
\boldsymbol {w} | _ {K} = \nabla (v | _ {K}), \quad \forall K \in \mathcal {M} _ {h},
$$

其中在单元边界上的值可以定义为任何有限值. 设 $\varphi \in C_0^\infty (\Omega)^d$ ，则

$$
\begin{array}{l} \int_ {\Omega} \boldsymbol {w} \cdot \varphi \mathrm {d} \boldsymbol {x} = \sum_ {K \in \mathcal {M} _ {h}} \int_ {K} \nabla v \cdot \varphi \mathrm {d} \boldsymbol {x} \\ = \sum_ {K \in \mathcal {M} _ {h}} \left(- \int_ {K} v \nabla \cdot \varphi d x + \int_ {\partial K} v \varphi \cdot n _ {K}\right) = - \int_ {\Omega} v \nabla \cdot \varphi d x, \\ \end{array}
$$

其中 $n_K$ 是 $\partial K$ 的单位外法向量.这里我们用到了 $v$ 的连续性及在任意公共面 $e = K\cap K^{\prime}$ 上 $n_K,n_{K^{\prime}}$ 反向.这说明 $\pmb{w}$ 是 $v$ 的弱梯度，因此 $v\in H^{1}(\varOmega)$

再证必要性. 设 $v \in H^{1}(\Omega)$ . 设 $e$ 是两个单元 $K_{1}$ 和 $K_{2}$ 的公共面, $\pmb{x}$ 是 $e$ 的一个内点. 则存在以 $\pmb{x}$ 为心的足够小的球 $B$ 使得 $B \subset K_{1} \cup K_{2}$ , 如图6.2所示. 记 $v_{i} = v|_{K_{i}}$ , $i = 1,2$ . 由Green公式, 对任意 $\varphi \in C_0^\infty(B)^d$ , 将 $\varphi$ 零延拓至 $\Omega$ 上有定义,

$$
\int_ {K _ {i}} \nabla v \cdot \varphi \mathrm {d} x = - \int_ {K _ {i}} v \nabla \cdot \varphi \mathrm {d} x + \int_ {\partial K _ {i}} v _ {i} (\varphi \cdot n _ {K _ {i}}), \quad i = 1, 2.
$$

既然 $v\in H^{1}(\varOmega)$ ，我们有

$$
\int_ {\Omega} \nabla v \cdot \varphi \mathrm {d} x = - \int_ {\Omega} v \nabla \cdot \varphi \mathrm {d} x.
$$

因此

$$
\int_ {e} \left(v _ {1} - v _ {2}\right) \phi \mathrm {d} s = 0, \quad \forall \phi \in C _ {0} ^ {\infty} (B).
$$

可推出 $v_{1}(\pmb {x}) = v_{2}(\pmb {x})$ .因此， $\mathcal{V}$ 在任何内部面上连续，从而 $v\in C(\varOmega)$ 证毕

![](images/b9fe70af0634e57fa3ccd8acc47d3bcca56829593dda09f5d3fb1031c9f9245d.jpg)  
图6.2 单元 $K_{1}$ 和 $K_{2}$ 公共边上一点 $\pmb{x}$ 及其邻域 $B$

所以，分片充分光滑的试探函数空间 $U_{h} \subset H^{1}(\Omega)$ ，需要其中的函数整体连续。而 $U_{h} \subset H^{2}(\Omega)$ 需要其中的函数整体连续可微。有限元函数的整体光滑性，往往通过自由度的选取可以得到保证。比如，上一节的分段线性函数空间 $U_{h}$ ，自由度选在了区间端点的函数值，自然就保证了有限元函数的整体连续性，从而（6.3）中的 $U_{h} \subset H^{1}(\Omega)$ 。

下面, 我们介绍几种具体的有限元和相应的有限元空间. 对自然数 $p \geqslant 0$ , 记 $P_{p}(\Omega)$ 为 $\Omega$ 上次数小于等于 $p$ 的多项式集合, 记 $Q_{p}(\Omega)$ 为 $\Omega$ 上按每个坐标分量的次数都小于等于 $p$ 的多项式集合. 例如, 当 $d = 2$ 时, $Q_{1}$ 为双线性函数空间. 显然, 对一维情形有 $Q_{p} = P_{p}$ .

# 6.2.2 一维高次元

为了提高有限元法的精度, 需要增加试探函数空间 $U_{h}$ 的维数. 这有两个途径, 一是加密网格剖分使单元最大直径 $h$ 变小, 节点参数 $\{u_{i}\}$ 增加. 二是增加分段多项式的次数, 这就是本节要介绍的高次元. 引进高次元是有限元法的重要技巧.

一次元是分段一次多项式，在每一单元 $K_{i} = [x_{i - 1}, x_{i}]$ 上含有两个待定系数，自由度是2,恰好由两个端点值决定.分段二次、三次及更高次多项式在每一单元上的自由度

增加了, 应当按哪种插值去确定它们呢? 一种是 Lagrange 型, 在单元内部增加插值节点; 另一种是 Hermite 型, 在节点引进高阶导数. 无论用哪一种插值, 都要求它们在整个区域上有一定的光滑度 (参见定理 6.4), 以保证双线性形式有意义.

$p$ 次Lagrange元

设 $p \geqslant 1$ . 先定义三元组 $(K, \mathcal{P}, \mathcal{N})$ . 取 $K$ 为小单元 $K_{i} = [x_{i-1}, x_{i}]$ , $\mathcal{P} = P_{p}(K)$ . 注意到一个次数小于等于 $p$ 一元多项式有 $p + 1$ 个自由度, 需要 $p + 1$ 个不同的点就可以确定它. 在 $K$ 上取 $p + 1$ 个不同的节点 (包括 $K$ 的两个端点)

$$
x _ {i - 1} = s _ {0} <   s _ {1} <   \dots <   s _ {p} = x _ {i}.
$$

取 $p + 1$ 个自由度为 $N_{j}(v) = v(s_{j}), j = 0,1,\dots ,p.$

显然, $K$ 和 $\mathcal{P}$ 分别满足有限元定义的 (i) 和 (ii). 由于 (iii) 与 $(\mathrm{iii})'$ 等价, 下面我们验证 $(\mathrm{iii})'$ . 设 $v \in \mathcal{P}$ 满足 $N_{j}(v) = v(s_{j}) = 0, j = 0,1,\dots,p$ , 即一个次数小于等于 $p$ 的多项式 $v$ 有 $p + 1$ 个不同的根, 所以只能 $v = 0$ , 从而 $\mathcal{N} = \{N_0,N_1,\dots,N_p\}$ 满足 $(\mathrm{iii})'$ . 这样我们就证明了 $(K,\mathcal{P},\mathcal{N})$ 是有限元.

下面推导节点基函数的公式. 注意到 $K$ 的位置和长度任意, 为了公式的简洁和便于程序实现, 像6.1节那样, 引入参考单元 $\hat{K} = [0,1]$ , 并通过仿射变换

$$
\xi = \frac {x - x _ {i - 1}}{h _ {i}} \tag {6.8}
$$

将单元 $K$ 变到 $\hat{K}$ 。我们先推导参考单元上的有限元节点基函数。记 $\hat{K}$ 上的节点为

$$
0 = \xi_ {0} <   \xi_ {1} <   \dots <   \xi_ {p} = 1,
$$

则可取 $s_j = x_{i - 1} + h_i\xi_j, j = 0,1,\dots ,p.$ 由节点基的定义6.2，即

$$
\text {求} \Phi_ {j} \in P _ {p} (\hat {K})   \text {使 得}   \Phi_ {j} (\xi_ {l}) = \delta_ {j l},   0 \leqslant j, l \leqslant p. \tag {6.9}
$$

这是 $p + 1$ 个Lagrange插值问题.所以 $\xi_{j},s_{j}$ 也称为插值节点.由Lagrange插值公式得参考单元上的节点基函数：

$$
\Phi_ {j} (\xi) = \prod_ {\substack {l = 0 \\ l \neq j}} ^ {p} \frac {\xi - \xi_ {l}}{\xi_ {j} - \xi_ {l}}, \quad j = 0, 1, \dots , p. \tag{6.10}
$$

插值节点最常用也是最简单的取法是取等距节点，即取 $\xi_{j} = \frac{j}{p}$ 此时

$$
\Phi_ {j} (\xi) = \prod_ {\substack {l = 0 \\ l \neq j}} ^ {p} \frac {p \xi - l}{j - l}, \quad j = 0, 1, \dots , p. \tag{6.11}
$$

简单计算即可给出线性、二次和三次等距节点Lagrange元在参考单元上的节点基函数，如表6.1.

![](images/16ff7c262707701541ffd80d5bc1e52151dfc1666f253f4249e7e4bd231ed2c8.jpg)

---

表 6.1 参考单元上等距节点 $p$ 次 Lagrange 元的节点基函数  

<table><tr><td>p</td><td>节点基函数</td></tr><tr><td>1</td><td>Φ0=1-ξ, Φ1=ξ</td></tr><tr><td>2</td><td>Φ0=(1-ξ)(1-2ξ), Φ1=4ξ(1-ξ), Φ2=ξ(2ξ-1)</td></tr><tr><td>3</td><td>Φ0=1/2(1-ξ)(1-3ξ)(2-3ξ), Φ1=9/2ξ(1-ξ)(2-3ξ), Φ2=9/2ξ(1-ξ)(3ξ-1), Φ3=1/2ξ(1-3ξ)(2-3ξ)</td></tr></table>

通过仿射变换 (6.8) 就可以将参考单元上的节点基函数变为单元 $K_{i}$ 上的节点基函数

$$
\phi_ {i, j} (x) = \Phi_ {j} \left(\frac {x - x _ {i - 1}}{h _ {i}}\right), \quad j = 0, 1, \dots , p. \tag {6.12}
$$

显然，任一有限元函数 $v \in P_p(K_i)$ 都可以表示为

$$
v (x) = \sum_ {j = 0} ^ {p} v \left(s _ {j}\right) \phi_ {i, j} (x), \quad x \in K _ {i}. \tag {6.13}
$$

给定 $\Omega$ 的网格剖分 $\mathcal{M}_h$ ，将每个单元上的有限元拼起来就可以构造试探函数空间。比如为求解6.1节的两点边值问题(6.2)，类似于(6.3)，我们可以如下定义 $p$ 次Lagrange有限元空间

$$
U _ {h} = \left\{v _ {h} \in C (\Omega): v _ {h} | _ {K _ {i}} \in P _ {p} (K _ {i}), i = 1, 2, \dots , n, v _ {h} (0) = 0 \right\}. \tag {6.14}
$$

由于每个单元的端点都是插值节点, 所以 $v_{h}$ 的连续性是很容易得到保证的, 只需对每个端点 $x_{i} = K_{i} \cap K_{i + 1}$ , 两个单元的公共自由度 $(v_{h}|_{K_{i}})(x_{i})$ 和 $(v_{h}|_{K_{i + 1}})(x_{i})$ 取成一样即可. 由定理6.4, 这样定义的试探函数空间 $U_{h} \subset H_{E}^{1}$ .

当然，将每个单元（小区间）上的节点和节点基适当“拼”起来就可以得到有限元空间 $U_{h}$ 的节点和节点基：

$$
\begin{array}{l} x _ {(i - 1) p + j} ^ {p} = x _ {i - 1} + h _ {i} \xi_ {j}, \quad i = 1, 2, \dots , n, j = 1, 2, \dots , p. \\ \phi_ {(i - 1) p + j} ^ {p} (x) = \left\{ \begin{array}{l l} \phi_ {i, j} (x), & x \in [ x _ {i - 1}, x _ {i} ], \\ 0, & \text {其 他 ,} \end{array} \quad i = 1, 2, \dots , n, j = 1, 2, \dots , p - 1, \right. \\ \phi_ {i p} ^ {p} (x) = \left\{ \begin{array}{l l} \phi_ {i, p} (x), & x \in [ x _ {i - 1}, x _ {i} ], \\ \phi_ {i + 1, 0} (x), & x \in [ x _ {i}, x _ {i + 1} ], \quad i = 1, 2, \dots , n - 1, \\ 0, & \text {其 他}, \end{array} \right. \\ \phi_ {n p} ^ {p} (x) = \left\{ \begin{array}{l l} \phi_ {n, p} (x), & x \in [ x _ {n - 1}, x _ {n} ], \\ 0, & \text {其 他}, \end{array} \right. \\ \end{array}
$$

图6.3给出了Lagrange二次元空间节点基函数的示意图

![](images/2c141b85d8acc5eb16475cfd88dd445cde4ef0bc1a069b47215a8af9a86192e9.jpg)

![](images/2d5aaed3c2333d3807d40482e1dd59180bdc441446996a80ed87e384f674c1d1.jpg)  
图6.3Lagrange二次元空间节点基函数曲线图

显然,

$$
\phi_ {l} ^ {p} \left(x _ {m} ^ {p}\right) = \delta_ {l m}, \quad l, m = 1, 2, \dots , n p,
$$

且任一 $v_{h}\in U_{h}$ 可以按节点基 $\{\phi_l^p,l = 1,2,\dots ,np\}$ 表示为

$$
v _ {h} (x) = \sum_ {l = 1} ^ {n p} v _ {h} \left(x _ {l} ^ {p}\right) \phi_ {l} ^ {p} (x).
$$

# 三次Hermite元

要想得到 $C^1 (\varOmega)$ 类试探函数空间，可采用分段三次多项式函数.先在每个单元上定义有限元.同样地，取 $K = K_{i} = [x_{i - 1},x_{i}]$ .再令 $\mathcal{P} = P_3(K)$ ．三次多项式有四个待定系数，需要4个自由度，我们分别取两个端点的函数值和一阶导数值.若记 $s_0 = x_{i - 1}$ $s_1 = x_i$ ，即取 $N_{00}(v) = v(s_0),N_{01}(v) = v'(s_0),N_{10}(v) = v(s_1),N_{11}(v) = v'(s_1)$ 作为自由度，令 $\mathcal{N} = \{N_{00},N_{01},N_{10},N_{11}\}$ ，则容易证明这样定义的 $(K,\mathcal{P},\mathcal{N})$ 是有限元.由于相应的有限元插值是Hermite型插值，所以称这个有限元是三次Hermite元.

下面考虑节点基函数. 同样地, 我们先在参考单元 $\hat{K} = [0,1]$ 上推导节点基函数公式. 即求 $\Phi_{j}^{(l)} \in P_{3}(\hat{K})$ 满足

$$
\Phi_ {j} ^ {(0)} (m) = \delta_ {j m}, \Phi_ {j} ^ {(0) \prime} (m) = 0; \quad \Phi_ {j} ^ {(1)} (m) = 0, \Phi_ {j} ^ {(1) \prime} (m) = \delta_ {j m}, \quad j, m = 0, 1. \tag {6.15}
$$

先求 $\Phi_0^{(0)}$ ，即求三次多项式使得

$$
\Phi_ {0} ^ {(0)} (0) = 1, \quad \Phi_ {0} ^ {(0) \prime} (0) = 0, \quad \Phi_ {0} ^ {(0)} (1) = \Phi_ {0} ^ {(0) \prime} (1) = 0.
$$

显然，1是二重根， $\Phi_0^{(0)}(\xi)$ 形如

$$
\Phi_ {0} ^ {(0)} (\xi) = (1 - \xi) ^ {2} (\alpha \xi + \beta).
$$

为使前两个条件满足，应取 $\alpha = 2,\beta = 1.$ 于是 $\varPhi_0^{(0)}(\xi)=(1-\xi)^2(2\xi+1)$ 再求三次式 $\varPhi_0^{(1)}$ 满足条件

$$
\Phi_ {0} ^ {(1)} (0) = 0, \quad \Phi_ {0} ^ {(1) \prime} (0) = 1, \quad \Phi_ {0} ^ {(1)} (1) = \Phi_ {0} ^ {(1) \prime} (1) = 0.
$$

显然， $\Phi_0^{(1)}(\xi)$ 形如

$$
\Phi_ {0} ^ {(1)} (\xi) = c \xi (1 - \xi) ^ {2}.
$$

选取常数 $c = 1$ 使得第二个条件成立即可. $\Phi_1^{(0)}$ 和 $\Phi_1^{(1)}$ 可类似求得. 总之, 我们可得到参考单元上的节点基函数:

$$
\begin{array}{l} \Phi_ {0} ^ {(0)} (\xi) = (1 - \xi) ^ {2} (2 \xi + 1), \quad \Phi_ {0} ^ {(1)} (\xi) = \xi (1 - \xi) ^ {2}, \tag {6.16} \\ \Phi_ {1} ^ {(0)} (\xi) = \xi^ {2} (3 - 2 \xi), \quad \Phi_ {1} ^ {(1)} (\xi) = \xi^ {2} (\xi - 1). \\ \end{array}
$$

利用上面的形函数及复合函数求导法则可得小单元 $K_{i}$ 上的节点基函数：

$$
\phi_ {i, j} ^ {(0)} (x) = \Phi_ {j} ^ {(0)} \left(\frac {x - x _ {i - 1}}{h _ {i}}\right), \quad \phi_ {i, j} ^ {(1)} (x) = h _ {i} \Phi_ {j} ^ {(1)} \left(\frac {x - x _ {i - 1}}{h _ {i}}\right), \quad j = 0, 1. \tag {6.17}
$$

显然，任一有限元函数 $v \in P_3(K_i)$ 都可以表示为

$$
v (x) = \sum_ {j = 0} ^ {1} \left(v \left(x _ {i - 1 + j}\right) \phi_ {i, j} ^ {(0)} (x) + v ^ {\prime} \left(x _ {i - 1 + j}\right) \phi_ {i, j} ^ {(1)} (x)\right), \quad \forall x \in K _ {i}. \tag {6.18}
$$

以上方法同样可用来构造四次元、五次元以及更高次元

同样, 给定 $\Omega$ 的网格剖分 $\mathcal{M}_h$ , 我们可以引入如下三次Hermite有限元空间:

$$
U _ {h} = \left\{v _ {h} \in C ^ {1} (\Omega): v _ {h} | _ {K _ {i}} \in P _ {3} (K _ {i}), i = 1, 2, \dots , n, v _ {h} (0) = 0 \right\}. \tag {6.19}
$$

由于自由度取成了每个单元的端点处的函数值和导数值, 所以 $v_{h}$ 的连续可微性是很容易得到保证的. 由定理6.4, 这样定义的试探函数空间 $U_{h} \subset H^{2}$ .

类似于线性有限元方法 (6.4), 将其中的分段线性有限元空间换为 (6.19) 中的 $U_h$ , 即得求解两点边值问题 (6.1) 的三次 Hermite 有限元方法: 求 $u_h \in U_h$ 使得

$$
a \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in U _ {h}. \tag {6.20}
$$

为了推导相应的有限元方程组，我们需要将 $U_{h}$ 中的有限元函数按节点基展开.将 $\phi_{j,l}^{(i)}$ 在单元 $K_{i}$ 外作零延拓，仍记为 $\phi_{j,l}^{(i)}$ .给定 $v_{h}$ 在节点 $x_{i}$ 处的函数值 $v_{i} = v_{h}(x_{i})$ 和导数值 $v_{i}^{\prime} = v_{h}^{\prime}(x_{i})$ ，由(6.17)，知 $v_{h}$ 可表示为

$$
\begin{array}{l} v _ {h} (x) = \sum_ {i = 1} ^ {n} \sum_ {j = 0} ^ {1} \left(v _ {i - 1 + j} \phi_ {i, j} ^ {(0)} (x) + v _ {i - 1 + j} ^ {\prime} \phi_ {i, j} ^ {(1)} (x)\right) \\ = \sum_ {i = 0} ^ {n} \left(v _ {i} \left(\phi_ {i, 1} ^ {(0)} (x) + \phi_ {i + 1, 0} ^ {(0)} (x)\right) + v _ {i} ^ {\prime} \left(\phi_ {i, 1} ^ {(1)} (x) + \phi_ {i + 1, 0} ^ {(1)} (x)\right)\right), \tag {6.21} \\ \end{array}
$$

这里规定 $\phi_{0,1}^{(l)} = \phi_{n + 1,0}^{(l)}\equiv 0,l = 1,2.$ 记

$$
\varphi_ {i} ^ {(0)} = \phi_ {i, 1} ^ {(0)} (x) + \phi_ {i + 1, 0} ^ {(0)} (x), \quad \varphi_ {i} ^ {(1)} = \phi_ {i, 1} ^ {(1)} (x) + \phi_ {i + 1, 0} ^ {(1)} (x). \tag {6.22}
$$

易知 $\varphi_{i}^{(l)}(l = 0,1)$ 连续可微，且满足

$$
\varphi_ {i} ^ {(0)} (x _ {j}) = \delta_ {i j}, \varphi_ {i} ^ {(0) \prime} (x _ {j}) = 0; \quad \varphi_ {i} ^ {(1)} (x _ {j}) = 0, \varphi_ {i} ^ {(1) \prime} (x _ {j}) = \delta_ {i j},
$$

所以 $\varphi_i^{(l)}(l = 0,1)$ 是有限元空间 $U_{h}$ 的节点基函数

图6.4给出了Hermite三次元空间节点基函数的示意图

![](images/ef95aabbfc6078a519dd7e49568cb5546047a21cdca5614eb75eeb6ed800307d.jpg)

![](images/0436595b538028bd2f1637bf90350bc7447410f7501537a7c7dc8e299f91ed7a.jpg)  
图6.4Hermite三次元空间节点基函数

有限元解可表示为

$$
u _ {h} = \sum_ {i} \left[ u _ {i} \varphi_ {i} ^ {(0)} (x) + u _ {i} ^ {\prime} \varphi_ {i} ^ {(1)} (x) \right], \tag {6.23}
$$

其中 $u_{i} = u_{h}(x_{i}),u_{i}^{\prime} = u_{h}^{\prime}(x_{i})$ .在(6.20）中取检验函数 $v_{h} = \varphi_{i}^{(l)}$ ，即得有限元方程组

$$
\sum_ {i = 0} ^ {n} \left[ u _ {i} a \left(\varphi_ {i} ^ {(0)}, \varphi_ {j} ^ {(l)}\right) + u _ {i} ^ {\prime} a \left(\varphi_ {i} ^ {(1)}, \varphi_ {j} ^ {(l)}\right) \right] = \left(f, \varphi_ {j} ^ {(l)}\right), \quad l = 0, 1, \quad j = 1, 2, \dots , n. \tag {6.24}
$$

左端点方程为

$$
\left\{ \begin{array}{l} u _ {0} = 0, \\ \sum_ {i = 0} ^ {n} \left[ u _ {i} a \left(\varphi_ {i} ^ {(0)}, \varphi_ {0} ^ {(1)}\right) + u _ {i} ^ {\prime} a \left(\varphi_ {i} ^ {(1)}, \varphi_ {0} ^ {(1)}\right) \right] = (f, \varphi_ {0} ^ {(1)}). \end{array} \right. \tag {6.25}
$$

同线性有限元方法一样，可以先计算单元刚度矩阵，再组装总刚度矩阵。当然，工作量要大一些。

方程 (6.24) (6.25) 有两组未知量: $\{u_i\}$ 和 $\{u_i'\}$ (称为广义坐标), 相当于固体力学中的位移和应力, 都是实际要计算的量. 假若用一次元, 求出 $u_i$ 后, 还要作一次微商运算 (左、右微商) 和加权平均, 才得到 $u_i'$ . 这样, 往往会影响应力的精度. 用高次元则可改善应力的计算.

高次元的另一个优点是收敛阶高. 例如三次元, 若精确解足够光滑, 则收敛阶可达到 $O(h^{3})$ (按 $H^{1}$ 度量), 这就可以适当放大步长.

另一方面也要看到，采用高次元要付出一些代价.除了增加计算积分的复杂性外，刚度矩阵的带宽也比一次元更大了.

此外，若边值问题的解本身不够光滑，用高次元就不能达到高精度的目的.所以，用

哪种类型的试探函数作有限元逼近，要根据问题的性质和机器条件决定.

# 6.2.3 解二维问题的矩形元与四边形元

从本节起，我们讨论二、三维椭圆边值问题的有限元解法.首先讨论一些常用单元的形状函数及构造方法，包括矩形元、四边形元、三角形元、四面体元和长方体元.本节讨论矩形元和四边形元

Lagrange 型矩形元公式

为简单起见, 假定区域 $G$ 可以分割成有限个矩形的和, 且每个小矩形的边和坐标轴平行. 任意两个矩形, 或者不相交, 或者有公共的边或公共的顶点. 我们把每一小矩形叫做单元, 称如此的分割为矩形剖分.

取定剖分后，我们着手构造Lagrange有限元. 取 $K$ 为任一矩形单元， $K = K_{ij} = [x_{i-1}, x_i] \times [y_{j-1}, y_j]$ . 给定 $p \geqslant 1$ , 取 $\mathcal{P} = Q_p(K)$ 为 $K$ 上的双 $p$ 次多项式空间. 共需要 $(p+1)^2$ 个自由度. 将矩形沿 $x, y$ 方向都作 $p$ 等分, 记 $h_{xi} = x_i - x_{i-1}$ , $h_{yj} = y_j - y_{j-1}$ 为网格步长. 两个方向的等分点记为 $s_l = x_{i-1} + lh_{xi}$ , $t_m = y_{j-1} + mh_{yj}$ ( $l, m = 0, 1, \dots, p$ ). 取插值节点为 $(s_l, t_m)$ (如图6.5), 则 $(p+1)^2$ 个自由度定义为 $N_{lm}(v) = v(s_l, t_m)$ , $l, m = 0, 1, \dots, p$ .

![](images/172a13ecf3295246b6d563a23bfe6a5eb95c70c908367deda365722517e734e7.jpg)

![](images/efb99e5e7c6e00586313acfae4d1c04e68b1fe435f3764062e61e1518c965be4.jpg)

![](images/6652a16374979f5a47de0679786a9b2d263d6690b7a3d73b007ec0fb54dc7d40.jpg)  
图6.5 Lagrange双 $p$ 次元插值节点分布， $p = 1,2,3$

下面验证有限元的定义. 只需验证 $(\mathrm{iii})'$ . 假设 $v \in Q_p(K)$ 且 $N_{lm}(v) = v(s_l, t_m) = 0$ , $l, m = 0, 1, \dots, p$ . 由 $v(s_l, y)$ 是 $y$ 的 $p$ 次式, 得 $v(s_l, y) \equiv 0$ , 故 $v(x, y)$ 能被 $x - s_l$ 整除. 同理 $v$ 也能被 $y - t_m$ 整除. 即 $v = c \prod_{l=0}^{p} \prod_{m=0}^{p} (x - s_l)(y - t_m)$ , 得 $v \equiv 0$ . 所以上面定义的三元组是一个有限元, 称为 Lagrange 双 $p$ 次元.

为了给出单元 $K_{ij}$ 上的节点基函数, 同一维情形一样, 我们首先推导参考单元 $\hat{K} = [0,1] \times [0,1]$ 上的节点基函数. 插值节点为 $(\xi_l, \xi_m), \xi_l = \frac{l}{p}, l, m = 0,1,\dots,p.$ 利用一维参考单元上的 Lagrange 形函数 $\Phi_j(\xi)$ (见 (6.11) 及表 6.1) 作乘积, 就可以得到 Lagrange 双 $p$ 次元在参考单元上的节点基函数:

$$
\Phi_ {l m} (\xi , \eta) = \Phi_ {l} (\xi) \Phi_ {m} (\eta), \quad l, m = 0, 1, \dots , p. \tag {6.26}
$$

由(6.9)，上面定义的 $\varPhi_{lm}$ 确实满足：

$$
\Phi_ {l m} \left(\xi_ {l ^ {\prime}}, \xi_ {m ^ {\prime}}\right) = \delta_ {l l ^ {\prime}} \delta_ {m m ^ {\prime}}, \quad 0 \leqslant l, l ^ {\prime}, m, m ^ {\prime} \leqslant p.
$$

一旦得到了参考单元上的节点基函数，通过仿射变换

$$
\xi = \frac {x - x _ {i - 1}}{h _ {x i}}, \quad \eta = \frac {y - y _ {j - 1}}{h _ {y j}}, \tag {6.27}
$$

就可以得到小单元 $K_{ij}$ 上的节点基函数：

$$
\phi_ {i j, l m} (x, y) = \Phi_ {l} \left(\frac {x - x _ {i - 1}}{h _ {x i}}\right) \Phi_ {m} \left(\frac {y - y _ {j - 1}}{h _ {y j}}\right). \tag {6.28}
$$

显然，任一有限元函数 $v \in Q_p(K_{ij})$ 都可以表示为

$$
v (x, y) = \sum_ {l, m = 0} ^ {p} v \left(s _ {l}, t _ {m}\right) \phi_ {i j, l m} (x, y), \quad (x, y) \in K _ {i j}. \tag {6.29}
$$

把每个单元上的有限元拼起来, 就可以得到矩形网上相应的双 $p$ 次Lagrange有限元空间 $U_{h} \subset H^{1}$ .

# Hermite 型矩形元公式

矩形单元上一种有代表性的Hermite型有限元为Bonger-Fox-Schmit元，简称BFS元.设 $K = K_{ij} = [x_{i - 1},x_i]\times [y_{j - 1},y_j]$ 为矩形剖分中的任一矩形单元，其四个顶点分别记为 $P_{1},P_{2},P_{3},P_{4}$ ，如图6.6所示．取 $\mathcal{P} = Q_3(K)$ 为 $K$ 上的双三次多项式空间．16个自由度取为矩形单元四个顶点的函数值、8个一阶偏导数值和4个混合二阶偏导数值，即

$$
N _ {P _ {k}} ^ {(1)} (v) = v \left(P _ {k}\right),
$$

$$
N _ {P _ {k}} ^ {(2)} (v) = \frac {\partial v}{\partial x} (P _ {k}),
$$

$$
N _ {P _ {k}} ^ {(3)} (v) = \frac {\partial v}{\partial y} (P _ {k}),
$$

$$
N _ {P _ {k}} ^ {(4)} (v) = \frac {\partial^ {2} v}{\partial x \partial y} (P _ {k}),
$$

其中 $k = 1,2,3,4.$ 令 $\mathcal{N} = \{N_{P_k}^{(l)}\}_{k,l = 1}^4$ 容易验证这样定义的 $(K,\mathcal{P},\mathcal{N})$ 是有限元

下面给出BFS元的节点基函数.事实上，BFS元可以看作是一维Hermite元的张量积形式，记

$$
i _ {1} = 0, i _ {2} = 1, i _ {3} = 0, i _ {4} = 1,
$$

$$
j _ {1} = 0, j _ {2} = 0, j _ {3} = 1, j _ {4} = 1,
$$

![](images/8e3f8f845bc491ec2b5621bdfd2ed2b0e44b1b27ae59d2257b7dec33b5acdafe.jpg)  
图6.6BFS元在单元 $K$ 上的自由度

则根据 (6.16) 式就可得到单元 $\hat{K}$ 上的节点基函数：

$$
\Phi_ {k} ^ {(l)} (\xi , \eta) = \Phi_ {i _ {k}} ^ {(i _ {l})} (\xi) \Phi_ {j _ {k}} ^ {(i _ {l})} (\eta), \quad k, l = 1, 2, 3, 4.
$$

利用这些基函数及复合函数求导法则可得单元 $K$ 上的节点基函数

$$
\phi_ {P _ {k}} ^ {(l)} = h _ {x i} ^ {i _ {l}} h _ {y j} ^ {j _ {l}} \Phi_ {k} ^ {(l)} \left(\frac {x - x _ {i - 1}}{h _ {x i}}, \frac {y - y _ {j - 1}}{h _ {y j}}\right), \quad k, l = 1, 2, 3, 4.
$$

设 $u \in Q_3(K)$ , 且 $u_{P_k}^{(1)}, u_{P_k}^{(2)}, u_{P_k}^{(3)}, u_{P_k}^{(4)}$ 分别表示函数 $u, \frac{\partial u}{\partial x}, \frac{\partial u}{\partial y}$ 和 $\frac{\partial^2 u}{\partial x \partial y}$ 在节点 $P_k$ 处的函数值, 则 $u$ 可唯一表示为

$$
u (x, y) = \sum_ {k, l = 1} ^ {4} u _ {P _ {k}} ^ {(l)} \phi_ {P _ {k}} ^ {(l)}, \quad (x, y) \in K. \tag {6.30}
$$

Bogner-Fox-Schmit元相应的有限元空间为

$$
U _ {h} = \left\{u _ {h} \in C ^ {1} (\Omega): u _ {h} | _ {K} \in Q _ {3} (K), \forall K \in \mathcal {T} _ {h} \right\}. \tag {6.31}
$$

显然试探函数空间 $U_{h}\subset H^{2}(\varOmega)$

四边形元

如果区域 $\Omega$ 的边界 $\Gamma$ 过于复杂，以致用折线逼近的几何误差太大，就需采取分段高次曲线逼近，这时将出现曲边单元.前面讨论了平行坐标轴的矩形剖分及矩形单元的形状函数，这类剖分对规则区域才是方便的；否则可采取任意四边形剖分，它和三角剖分一样有很大的灵活性.对于上述曲边单元、四边形单元，如何构造单元形状函数？本节我们以一般四边形元为例给出构造方法

前面构造单元形状函数时, 是用一个仿射变换 (6.27), 把任意矩形单元 $K$ 变到 $\xi \eta$ 平面上的“参考元” $\hat{K}$ , $\hat{K}$ 是单位正方形 $0 \leqslant \xi \leqslant 1, 0 \leqslant \eta \leqslant 1$ (参看图 6.7(a)), 其几何形状非常简单, 易于构造单元形状函数. 其实 $K$ 也是参考元 $\hat{K}$ 在仿射变换下的映像. 不同的仿射变换, 便得到不同的单元及其形状函数. 所以, $\Omega$ 的任一网格剖分及试探函数, 可以看作是参考元 $\hat{K}$ 及其形状函数在一族仿射变换下得到的.

例6.1再看一下 $p = 1$ 时的Lagrange矩形元.设参考有限元取为 $(\hat{K},Q_1(\hat{K}),\hat{\mathcal{N}})$ 即参考单元上的Lagrange双线性元，其中

$$
Q _ {1} (\hat {K}) = \left\{\left(a _ {0} + a _ {1} \xi\right) \left(b _ {0} + b _ {1} \eta\right): a _ {0}, a _ {1}, b _ {0}, b _ {1} \in \mathbb {R} \right\},
$$

自由度为参考单元四个顶点的函数值

$$
\hat {N} _ {1} (\hat {v}) = \hat {v} (0, 0), \hat {N} _ {2} (\hat {v}) = \hat {v} (1, 0), \hat {N} _ {3} (\hat {v}) = \hat {v} (0, 1), \hat {N} _ {4} (\hat {v}) = \hat {v} (1, 1).
$$

考虑矩形单元 $K = \square P_{1}P_{2}P_{3}P_{4}$ ，其中顶点 $P_{i}$ 为 $(x_{i},y_{i}),i = 1,2,3,4,$ 满足 $x_{3} = x_{1},x_{4} =$ $x_{2},y_{2} = y_{1},y_{4} = y_{3}.$ 通过仿射变换(参见(6.27))

$\left\{ \begin{array}{ll}x = x_1 + (x_2 - x_1)\xi ,\\ y = y_1 + (y_3 - y_1)\eta , \end{array} \right.$ 记为 $(x,y) = F_K(\xi ,\eta),$ (6.32)

可以将 $\hat{K}$ 变为 $K$ ，即 $K = F_{K}(\hat{K})$ ，将 $Q_{1}(\hat{K})$ 变为 $Q_{1}(K)$ ，即 $Q_{1}(K) = \{v:v\circ F_{K}\in$ $Q_{1}(\hat{K})\}$ ，四个自由度变为 $N_{i}(v) = v(P_{i}),i = 1,2,3,4$ (参看图6.7(b)).

![](images/2f3577bb335a2baa38fc900ecd5ed74abfa83c2d37598eab79e50ec9a0e2da87.jpg)  
(a) 参考单元为单位正方形

![](images/00813be643b7067752719742c01392add115e6e4946493b1eb8a7163c20d2530.jpg)  
(b) $\widehat{K}$ 变为矩形 $K$   
图6.7 参考单元

现在我们面临的是一般单元 $K$ (四边形单元、曲边单元), 只限于仿射变换是不够的. 因此需要考虑更一般的可逆连续变换:

$$
(x, y) = F _ {K} (\xi , \eta) = (x (\xi , \eta), y (\xi , \eta)), \tag {6.33}
$$

其中 $(\xi, \eta) \in \hat{K}, (x, y) \in K.$ 变换 (6.33) 和单元 $K$ 有关, 不同的 $K$ 对应不同的变换. 也可反过来说, 不同的变换对应不同的单元 $K$ . 变换 (6.33) 应满足下列要求: 第一, 具有必要的光滑性. 通常取它为 $\xi, \eta$ 的多项式, 所以光滑性条件恒满足. 第二, $F_K$ 应是 $\hat{K}$ 到 $K$ 的一对一的变换, 就是说, $F_K$ 的 Jacobi 行列式

$$
J (\xi , \eta) = \frac {\partial (x , y)}{\partial (\xi , \eta)} \neq 0, \quad (\xi , \eta) \in \hat {K}. \tag {6.34}
$$

应指出的是， $\hat{K}$ 上的形状函数 $p_{\hat{K}}(\xi ,\eta)$ 虽然是多项式，但是通过变换 $F_{K}$ 的逆变换

消去 $\xi, \eta$ 后就不一定是多项式了，它可能是有理函数，也可能是无理函数。好在我们并不需要消去 $\xi, \eta$ 得到以 $x, y$ 表示的形状函数，因为形成有限元方程时遇到的积分可通过变换（6.33）化为 $\hat{K}$ 上对 $\xi, \eta$ 的积分。

构造变换 (6.33) 的方法很多, 应用中最重要的一种是取它和形状函数具同样形式, 这就是所谓“等参变换”. 当然, 前面的仿射变换也是等参变换. 使用等参变换也是有限元法的一个技巧.

例6.2（任意四边形单元）如例6.1, 参考有限元取为 $(\hat{K}, Q_1(\hat{K}), \hat{\mathcal{N}})$ . 如图6.8, $K$ 是 $xy$ 平面上任一四边形, 顶点是 $P_i(x_i, y_i), i = 1, 2, 3, 4$ . 我们按照例6.1的流程定义 $K$ 上的有限元. 首先, 易知, 把 $\hat{K}$ 变到 $K$ 的变换 (6.33) 可定义为

$$
\left[ \begin{array}{c} x \\ y \end{array} \right] = \varPhi_ {0} (\xi) \varPhi_ {0} (\eta) \left[ \begin{array}{c} x _ {1} \\ y _ {1} \end{array} \right] + \varPhi_ {1} (\xi) \varPhi_ {0} (\eta) \left[ \begin{array}{c} x _ {2} \\ y _ {2} \end{array} \right] + \varPhi_ {0} (\xi) \varPhi_ {1} (\eta) \left[ \begin{array}{c} x _ {3} \\ y _ {3} \end{array} \right] + \varPhi_ {1} (\xi) \varPhi_ {1} (\eta) \left[ \begin{array}{c} x _ {4} \\ y _ {4} \end{array} \right],
$$

其中 $\varPhi_0(s)=1-s,\varPhi_1(s)=s.$ 或改写为

$$
\left\{ \begin{array}{l} x (\xi , \eta) = x _ {1} + \left(x _ {2} - x _ {1}\right) \xi + \left(x _ {3} - x _ {1}\right) \eta + \left(x _ {4} - x _ {3} - x _ {2} + x _ {1}\right) \xi \eta , \\ y (\xi , \eta) = y _ {1} + \left(y _ {2} - y _ {1}\right) \xi + \left(y _ {3} - y _ {1}\right) \eta + \left(y _ {4} - y _ {3} - y _ {2} + y _ {1}\right) \xi \eta . \end{array} \right. \tag {6.35}
$$

与 $\hat{K}$ 上的形状函数一样都是双线性函数，所以这个变换是一种等参变换，且当四边形 $K$ 的对角线中点不重合时，不是仿射变换。显然，变换 (6.35) 把 $\hat{K}$ 的每条边仿射变换到 $K$ 的对应边，把 $\hat{K}$ 的内部变到 $K$ 的内部。

为了检验变换 (6.35) 是一对一的, 计算 Jacobi 行列式

$$
J = \left| \begin{array}{c c} x _ {2} - x _ {1} + A \eta & x _ {3} - x _ {1} + A \xi \\ y _ {2} - y _ {1} + B \eta & y _ {3} - y _ {1} + B \xi \end{array} \right|,
$$

其中 $A = x_{4} - x_{3} - x_{2} + x_{1}, B = y_{4} - y_{3} - y_{2} + y_{1}$ . 展开后, $J$ 是 $\xi, \eta$ 的双线性函数. 故只需检验 $J$ 在四个顶点的值有相同符号即可. 于 $\xi = 0, \eta = 0$ ,

$$
J (0, 0) = \left(x _ {2} - x _ {1}\right) \left(y _ {3} - y _ {1}\right) - \left(y _ {2} - y _ {1}\right) \left(x _ {3} - x _ {1}\right) = l l ^ {\prime} \sin \theta ,
$$

长度 $l, l'$ 和角 $\theta$ 如图6.8所示. 当 $0 < \theta < \pi$ 时, $J(0,0) > 0$ . 同样, 其他三个内角小于 $\pi$ 时, $J$ 在 $(1,0), (0,1), (1,1)$ 也大于0. 总之, $J$ 于 $\hat{K}$ 恒大于0的充要条件是 $K$ 为凸四边形.

![](images/eb83896fb1db90aaa5f4caad3deabff40f93fe9bb9de6620ba37bcf75c8b477d.jpg)  
图6.8 等参变换

有了等参变换之后, 即可如下定义一般四边形 $K$ 上的有限元 $(K, \mathcal{P}, \mathcal{N})$ :

$$
\mathcal {P} = \left\{v: v \circ F _ {K} \in Q _ {1} (\hat {K}) \right\}, \quad N _ {i} (v) = v (P _ {i}), i = 1, 2, 3, 4.
$$

称为映射的双线性函数空间. 注意 (6.35) 的逆变换一般是无理函数, 所以 $\hat{K}$ 上的形状函数尽管是简单的 (双线性函数), 但对应到 $K$ 上则是 $x, y$ 的无理函数. 好在我们不必通过 $x, y$ 把 $K$ 上的形状函数以显式表示出来 (参看 6.3 节).

# 6.2.4 三角形元

对曲边区域 $\Omega$ ，一般采用三角网近似。不妨设 $\Omega$ 是多边形域（否则可用多边形域逼近它）。将 $\Omega$ 分割成有限个三角形之和，使不同三角形无重叠的内部，且任一三角形的顶点不属于其他三角形边的内部。这样就把 $\Omega$ 分割成三角形网，称为 $\Omega$ 的三角剖分。每个三角形称为单元，它的顶点称为节点。属同一单元的二顶点称为相邻节点，有公共边的两个三角形称为相邻单元。

由于三角剖分可构造非均匀网格，并且能较好地逼近具复杂边界的区域，所以在二维问题中，三角形元是应用最广的单元.

面积坐标及有关公式

设 $\triangle (i,j,k)$ 是以 $i,j,k$ 为顶点的任意三角形单元，面积为 $S$ 。我们约定 $i,j,k$ 的次序按逆时针方向排列。在 $\triangle (i,j,k)$ 内任取一点 $P$ ，坐标为 $(x,y)$ 。连接 $P$ 点与三个顶点，将 $\triangle (i,j,k)$ 分成三个三角形（参见图6.9）： $\triangle (i,j,P)$ ， $\triangle (j,k,P)$ ， $\triangle (k,i,P)$ ，其面积分别为 $S_{k}, S_{i}, S_{j}$ 。显然 $S_{i} + S_{j} + S_{k} = S$ 。令

$$
L _ {i} = \frac {S _ {i}}{S}, \quad L _ {j} = \frac {S _ {j}}{S}, \quad L _ {k} = \frac {S _ {k}}{S}, \tag {6.36}
$$

![](images/41f123ef4e560d406124e0aa58bc6cf80c447d78ac201854e983a3042b10ea0d.jpg)  
图6.9 面积坐标

则 $L_{i}, L_{j}, L_{k} \geqslant 0, L_{i} + L_{j} + L_{k} = 1$ . 给定一点 $P$ , 唯一确定如此的一组数 $(L_{i}, L_{j}, L_{k})$ . 反之, 任给一组 $(L_{i}, L_{j}, L_{k}), L_{i}, L_{j}, L_{k} \geqslant 0, L_{i} + L_{j} + L_{k} = 1$ , 按关系式 (6.36) 也唯一确定一点 $P$ . 所以同一点 $P$ , 既可用直角坐标 $(x, y)$ 表示, 也可用 $(L_{i}, L_{j}, L_{k})$ 表示. 我们称 $(L_{i}, L_{j}, L_{k})$ 为点 $P$ 的面积坐标. 因为三角形的面积与参考坐标系无关, 所以面积坐标也

与坐标系无关, 这是采用面积坐标的优点. 我们知道

$$
2 S = \left| \begin{array}{c c c} 1 & x _ {i} & y _ {i} \\ 1 & x _ {j} & y _ {j} \\ 1 & x _ {k} & y _ {k} \end{array} \right|, 2 S _ {i} = \left| \begin{array}{c c c} 1 & x & y \\ 1 & x _ {j} & y _ {j} \\ 1 & x _ {k} & y _ {k} \end{array} \right|, 2 S _ {j} = \left| \begin{array}{c c c} 1 & x _ {i} & y _ {i} \\ 1 & x & y \\ 1 & x _ {k} & y _ {k} \end{array} \right|, 2 S _ {k} = \left| \begin{array}{c c c} 1 & x _ {i} & y _ {i} \\ 1 & x _ {j} & y _ {j} \\ 1 & x & y \end{array} \right|.
$$

由此可建立面积坐标与直角坐标之间的下列转换关系：

$$
\begin{array}{l} \left\{ \begin{array}{l} L _ {i} = \frac {1}{2 S} \left[ \left(x _ {j} y _ {k} - x _ {k} y _ {j}\right) + \left(y _ {j} - y _ {k}\right) x + \left(x _ {k} - x _ {j}\right) y \right], \\ L _ {j} = \frac {1}{2 S} \left[ \left(x _ {k} y _ {i} - x _ {i} y _ {k}\right) + \left(y _ {k} - y _ {i}\right) x + \left(x _ {i} - x _ {k}\right) y \right], \\ L _ {k} = \frac {1}{2 S} \left[ \left(x _ {i} y _ {j} - x _ {j} y _ {i}\right) + \left(y _ {i} - y _ {j}\right) x + \left(x _ {j} - x _ {i}\right) y \right]. \end{array} \right. (6.37) \\ \left\{ \begin{array}{l} x = x _ {i} L _ {i} + x _ {j} L _ {j} + x _ {k} L _ {k}, \\ y = y _ {i} L _ {i} + y _ {j} L _ {j} + y _ {k} L _ {k}. \end{array} \right. (6.38) \\ \end{array}
$$

在推导后一关系式时，利用了等式 $L_{i} + L_{j} + L_{k} = 1$

由连锁规则不难看出

$$
\left\{ \begin{array}{l} \frac {\partial}{\partial x} = \frac {1}{2 S} \left[ \left(y _ {j} - y _ {k}\right) \frac {\partial}{\partial L _ {i}} + \left(y _ {k} - y _ {i}\right) \frac {\partial}{\partial L _ {j}} + \left(y _ {i} - y _ {j}\right) \frac {\partial}{\partial L _ {k}} \right], \\ \frac {\partial}{\partial y} = \frac {1}{2 S} \left[ \left(x _ {k} - x _ {j}\right) \frac {\partial}{\partial L _ {i}} + \left(x _ {i} - x _ {k}\right) \frac {\partial}{\partial L _ {j}} + \left(x _ {j} - x _ {i}\right) \frac {\partial}{\partial L _ {k}} \right]. \end{array} \right. \tag {6.39}
$$

利用 $L_{i} = 1 - L_{j} - L_{k}$ ，消去(6.38)右端的 $L_{i}$ ，则得到由 $L_{j}L_{k}$ 平面到 $xy$ 平面的仿射变换，其逆变换把 $\triangle (i,j,k)$ 变到 $L_{j}L_{k}$ 平面以 $(0,0),(0,1),(1,0)$ 为顶点的直角三角形，如图6.10，

$$
\begin{array}{l} \left(x _ {i}, y _ {i}\right) \leftrightarrow (0, 0), \\ (x _ {j}, y _ {j}) \leftrightarrow (1, 0), \\ (x _ {k}, y _ {k}) \leftrightarrow (0, 1). \\ \end{array}
$$

![](images/b5959520a412483b8554fc3a7e9d27d6762f0c500fe6d709712f77cdec31d258.jpg)  
图6.10 参考三角形单元

利用重积分变量替换公式，不难得出下列积分公式

$$
\iint_ {\triangle (i, j, k)} L _ {i} ^ {p} L _ {j} ^ {q} L _ {k} ^ {r} \mathrm {d} x \mathrm {d} y = 2 S \frac {p ! q ! r !}{p + q + r + 2}, \tag {6.40}
$$

其中 $p, q, r$ 是任意非负整数

Lagrange 型公式

取 $K$ 为任一三角形单元, 比如为了记号简单, $K = \triangle (1,2,3)$ . 给定 $p\geqslant 1$ , 令 $\mathcal{P} = P_p(K)$ . 易知 $\dim \mathcal{P} = n = \frac{1}{2} p(p + 1)$ . 需要 $n$ 个自由度. 将三角形 $K$ 的每条边 $p$ 等分, 连接任二相邻边的对应等分点作平行于另一条边的线段, 正好交于 $n$ 个点. 我们就取这 $n$ 个交点为插值节点 (图6.11画出了线性元和二次元的节点), 记为 $A_{1}, A_{2}, \dots, A_{n}$ , 则 $n$ 个自由度取为 $\mathcal{N} = \{N_i(v) = v(A_i), i = 1,2,\dots,n\}$ . 这样定义的有限元称为Lagrange $p$ 次元.

![](images/3166b4d627f00f8b8bdcb2b5a1734f08e46853ced0591823d1441421d76f21c1.jpg)  
(a) 三个顶点为插值节点

![](images/8c50552e4b5f8cdd3cf9833477a97c9d86b0c8bce2c6da30787f71dae1298cd6.jpg)  
(b) 三个顶点及三边中点为插值节点  
图6.11 线性元和二次元的插值节点

当 $p = 1$ 时，由不共线三点决定一个平面，知Lagrange线性元满足有限元的定义显然，三个面积坐标函数满足：

$$
L _ {i} (A _ {j}) = \delta_ {i j}, \quad i, j = 1, 2, 3.
$$

它们正好组成了线性元的节点基. 对任意 $v \in P_1(K)$ 有:

$$
v = v _ {1} L _ {1} + v _ {2} L _ {2} + v _ {3} L _ {3}, \quad v _ {i} = v (A _ {i}), i = 1, 2, 3. \tag {6.41}
$$

对 $p = 2$ ，我们验证条件 $(\mathrm{iii})'$ .假设 $v\in P_2(K),v(A_i) = 0,i = 1,2,\dots ,6.$ 由Cartesian坐标和面积坐标的关系，知 $v$ 是 $L_{1},L_{2}$ 的二次多项式，记为 $v(L_{1},L_{2})$ ，由 $v(A_{2}) = v(A_{3}) = v(A_{4}) = 0$ 得 $v(0,1) = v(0,0) = v\left(0,\frac{1}{2}\right) = 0.$ 再注意到 $v(0,L_2)$ 是 $L_{2}$ 的二次多项式，可得 $v(0,L_2)\equiv 0$ ，故 $L_{1}|v$ ，即 $L_{1}$ 整除 $v$ .同理， $L_{2}|v,L_{3}|v$ 即 $v = L_{1}L_{2}L_{3}w,w$ 是多项式.因为 $v$ 是二次式，所以 $w = 0$ ，从而 $v = 0$ .故Lagrange二次元满足有限元的定义.利用待定系数法及面积坐标，容易推出Lagrange二次元在三个顶点和三边中点的节点基函数可分别表示为

$$
L _ {i} (2 L _ {i} - 1), i = 1, 2, 3; \quad 4 L _ {i} L _ {j}, 1 \leqslant i <   j \leqslant 3.
$$

对任意 $v \in P_2(K)$ 有：

$$
v = \sum_ {i = 1} ^ {3} \left[ L _ {i} (2 L _ {i} - 1) v _ {i} + 4 L _ {i + 1} L _ {i + 2} v _ {3 + i} \right], \quad v _ {i} = v \left(A _ {i}\right), i = 1, 2, \dots , 6, \tag {6.42}
$$

其中 $L_{4} = L_{1},L_{5} = L_{2}$

还可以构造三次及高次的Lagrange型公式，但常用的是一次及二次插值公式

容易证明，由 $m$ 次Lagrange型插值公式生成的试探函数属于 $C$ ，但不属于 $C^1$ ，因此只能有 $U_{h}\subset H^{1}$

Hermite型公式

取 $K = \triangle (1,2,3)$ 考虑三次Hermite元，令 $\mathcal{P} = P_3(K)$ ： $\mathcal{N}$ 需要10个自由度，选取方法如图6.12所示，其中“0”表示重心

![](images/b699eeaa00a9f2cf336bbd62d0fac08ed33c1819e35c92f11c45f1e4a7ec0727.jpg)  
图6.12 插值节点

下面验证有限元的定义, 只需考虑条件 (iii)'. 记 $K$ 的三个顶点为 $A_{1}, A_{2}, A_{3}$ , 重心为 $A_{0}$ . 设 $v \in P_{3}(K)$ , $v(A_{i}) = v_{x}(A_{i}) = v_{y}(A_{i}) = v(A_{0}) = 0$ , $i = 1, 2, 3$ . 与 Lagrange 二次元的推导类似, 可以证明 $L_{i}|v, i = 1, 2, 3$ , 即 $v = L_{1}L_{2}L_{3}w$ , $w$ 是常数. 又因为 $v(A_{0}) = 0$ , 所以 $w = 0$ , 从而 $v = 0$ . 故三次 Hermite 元满足有限元的定义.

利用待定系数法及面积坐标, 可以推导出任意 $v \in P_3(K)$ 可表示为

$$
v = \alpha_ {0} ^ {(3)} v (A _ {0}) + \sum_ {i = 1} ^ {3} \left(\alpha_ {i} ^ {(3)} v (A _ {i}) + \beta_ {i} ^ {(3)} v _ {x} (A _ {i}) + \gamma_ {i} ^ {(3)} v _ {y} (A _ {i})\right), \tag {6.43}
$$

其中节点基函数：

$$
\begin{array}{l} \left\{ \begin{array}{l} \alpha_ {0} ^ {(3)} (x, y) = 2 7 L _ {1} L _ {2} L _ {3}, \\ \alpha_ {i} ^ {(3)} (x, y) = L _ {i} ^ {3} + 3 L _ {i} ^ {2} \left(L _ {j} + L _ {k}\right) - 7 L _ {i} L _ {j} L _ {k}, \end{array} \right. \tag {6.44} \\ \left\{ \begin{array}{l} \beta_ {i} ^ {(3)} (x, y) = \left(x _ {j} - x _ {i}\right) \left(L _ {i} ^ {2} L _ {j} - L _ {i} L _ {j} L _ {k}\right) + \left(x _ {k} - x _ {i}\right) \left(L _ {i} ^ {2} L _ {k} - L _ {i} L _ {j} L _ {k}\right), \\ \gamma_ {i} ^ {(3)} (x, y) = \left(y _ {j} - y _ {i}\right) \left(L _ {i} ^ {2} L _ {j} - L _ {i} L _ {j} L _ {k}\right) + \left(y _ {k} - y _ {i}\right) \left(L _ {i} ^ {2} L _ {k} - L _ {i} L _ {j} L _ {k}\right), \end{array} \right. \\ \end{array}
$$

其中 $j = i + 1, k = i + 2$

给定了一组广义坐标后（包括三角单元顶点的函数值，两个一阶偏导数及单元重心的函数值），由公式 (6.43)（让 $\triangle(1,2,3)$ 取遍一切单元）确定出整个剖分上的试探函数

$u_{h}(x,y)$ (下标 $h$ 表示一切单元的最大直径). 一切可能的试探函数构成试探函数空间 $U_{h}$ . 设 $K_{1}, K_{2}$ 是两相邻单元, $l$ 为其公共边, $u_{h}^{1}$ 和 $u_{h}^{2}$ 分别表示 $u_{h}$ 从 $K_{1}$ 和 $K_{2}$ 延拓到 $l$ 上的一元三次多项式. 因为 $u_{h}^{1}$ 和 $u_{h}^{2}$ 在 $l$ 的两端取相同函数值及沿 $l$ 方向的导数值 (它们由偏导数 $u_{x}, u_{y}$ 唯一确定), 故 $u_{h}^{1}$ 和 $u_{h}^{2}$ 在 $l$ 上恒等, 这说明 $u_{h} \in C$ , 因而 $U_{h} \subset H^{1}$ . 其次, 考察它们的一阶导数. 显然 $u_{h}^{1}$ 和 $u_{h}^{2}$ 沿 $l$ 方向的导数恒等, 但沿 $l$ 的法向导数不一定相等. 因为法向导数是 $l$ 上的一元二次多项式, 二者仅在端点相等, 不能完全确定该二次多项式. 因此 $u_{h}$ 一般不属于 $C^{1}$ , 从而也不能要求 $U_{h}$ 属于 $H^{2}$ .

应当指出，重心点的方程只包含重心值及它所属单元三个顶点的广义坐标，所以求解有限元方程时，可先消去重心处的广义坐标，降低方程组的阶，然后再解这个阶数较低的方程组.

# 6.2.5 三维有限元

本节中, 我们简要介绍三维空间中最简单的四面体线性元与长方体三线性元.

四面体线性元

设 $K = \triangle (1,2,3,4)$ 是以 $A_{1},A_{2},A_{3},A_{4}$ 为顶点的任意四面体单元.取 $\mathcal{P} = P_1(K)$ 令 $\mathcal{N} = \{N_i(v) = v(A_i),i = 1,2,3,4\}$ 称为 $K$ 上的Lagrange线性元

为了描述节点基函数，类似于二维三角形，可以对三维四面体引入体积坐标。在 $\triangle (1,2,3,4)$ 内任取一点 $P$ ，坐标为 $(x,y,z)$ 。连接 $P$ 点与四个顶点，将 $\triangle (1,2,3,4)$ 分成四个四面体： $\triangle (2,4,3,P)$ ， $\triangle (3,4,1,P)$ ， $\triangle (4,2,1,P)$ ， $\triangle (1,2,3,P)$ ，其体积分别为 $V_{1}, V_{2}, V_{3}, V_{4}$ 。显然 $V_{1} + V_{2} + V_{3} + V_{4} = V$ 。令

$$
L _ {1} = \frac {V _ {1}}{V}, \quad L _ {2} = \frac {V _ {2}}{V}, \quad L _ {3} = \frac {V _ {3}}{V}, \quad L _ {4} = \frac {V _ {4}}{V}, \tag {6.45}
$$

则 $L_{1}, L_{2}, L_{3}, L_{4} \geqslant 0, L_{1} + L_{2} + L_{3} + L_{4} = 1.$ 给定一点 $P$ ，唯一确定如此的一组数 $(L_{1}, L_{2}, L_{3}, L_{4})$ 。反之，任给一组 $(L_{1}, L_{2}, L_{3}, L_{4}), L_{1}, L_{2}, L_{3}, L_{4} \geqslant 0, L_{1} + L_{2} + L_{3} + L_{4} = 1.$ 按关系式 (6.45) 也唯一确定一点 $P$ 。所以同一点 $P$ ，既可用直角坐标 $(x, y, z)$ 表示，也可用 $(L_{1}, L_{2}, L_{3}, L_{4})$ 表示。我们称 $(L_{1}, L_{2}, L_{3}, L_{4})$ 为点 $P$ 的体积坐标。注意到在单元重心处体积坐标每个分量都相等，二维时面积坐标也满足同样性质，所以面积坐标和体积坐标也统称为重心坐标。

不妨设四面体四个顶点的排序满足右手螺旋规则, 则有

$$
V = \frac {1}{6} \left| \begin{array}{c c c c} 1 & x _ {1} & y _ {1} & z _ {1} \\ 1 & x _ {2} & y _ {2} & z _ {2} \\ 1 & x _ {3} & y _ {3} & z _ {3} \\ 1 & x _ {4} & y _ {4} & z _ {4} \end{array} \right|,
$$

由此可建立体积坐标与直角坐标之间的转换关系：

$$
\begin{array}{l} \left\{ \begin{array}{c c c c} L _ {1} = \frac {1}{6 V} \left| \begin{array}{c c c c} 1 & x & y & z \\ 1 & x _ {2} & y _ {2} & z _ {2} \\ 1 & x _ {3} & y _ {3} & z _ {3} \\ 1 & x _ {4} & y _ {4} & z _ {4} \end{array} \right|, & L _ {2} = \frac {1}{6 V} \left| \begin{array}{c c c c} 1 & x _ {1} & y _ {1} & z _ {1} \\ 1 & x & y & z \\ 1 & x _ {3} & y _ {3} & z _ {3} \\ 1 & x _ {4} & y _ {4} & z _ {4} \end{array} \right|, \\ L _ {3} = \frac {1}{6 V} \left| \begin{array}{c c c c} 1 & x _ {1} & y _ {1} & z _ {1} \\ 1 & x _ {2} & y _ {2} & z _ {2} \\ 1 & x & y & z \\ 1 & x _ {4} & y _ {4} & z _ {4} \end{array} \right|, & L _ {4} = \frac {1}{6 V} \left| \begin{array}{c c c c} 1 & x _ {1} & y _ {1} & z _ {1} \\ 1 & x _ {2} & y _ {2} & z _ {2} \\ 1 & x _ {3} & y _ {3} & z _ {3} \\ 1 & x & y & z \end{array} \right|. \end{array} \right. \\ \left\{ \begin{array}{l} x = x _ {1} L _ {1} + x _ {2} L _ {2} + x _ {3} L _ {3} + x _ {4} L _ {4}, \\ y = y _ {1} L _ {1} + y _ {2} L _ {2} + y _ {3} L _ {3} + y _ {4} L _ {4}, \\ z = z _ {1} L _ {1} + z _ {2} L _ {2} + z _ {3} L _ {3} + z _ {4} L _ {4}. \end{array} \right. \tag {6.46} \\ \end{array}
$$

显然，三维Lagrange线性元的节点基函数就是 $\{L_1,L_2,L_3,L_4\}$ ，且任一函数 $v\in$ $P_{1}(K)$ 可表示为

$$
v = v \left(A _ {1}\right) L _ {1} + v \left(A _ {2}\right) L _ {2} + v \left(A _ {3}\right) L _ {3} + v \left(A _ {4}\right) L _ {4}. \tag {6.47}
$$

易见，由此插值公式生成的试探函数整体属于 $C$ ，因此相应的线性有限元空间 $U_{h} \subset H^{1}$ 。

长方体三线性元

取 $K = K_{ijk} = [x_{i - 1},x_i]\times [y_{j - 1},y_j]\times [z_{k - 1},z_k]$ 为任一长方体单元， $\mathcal{P} = Q_1(K)$ 记 $K$ 的8个顶点为 $A_{1},A_{2},\dots ,A_{8}$ ，自由度取为 $N_{l}(v) = v(A_{l}),l = 1,2,\dots ,8.$ 显然此三元组满足有限元的定义，称为Lagrange三线性元.

类似于二维Lagrange矩形元，可以利用一维有限元节点基函数的乘积来定义三维长方体上有限元节点基函数.记 $h_{xi} = x_i - x_{i - 1},h_{yj} = y_j - y_{j - 1},h_{zk} = z_k - z_{k - 1}$ 为网格步长.通过仿射变换

$$
\xi = \frac {x - x _ {i - 1}}{h _ {x i}}, \quad \eta = \frac {y - y _ {j - 1}}{h _ {y j}}, \quad \zeta = \frac {z - z _ {k - 1}}{h _ {z k}}, \tag {6.48}
$$

总可将 $K$ 变成参考单元 $\hat{K} = [0,1]\times [0,1]\times [0,1]$ .如果在 $\hat{K}$ 上造出了节点基函数，再通过变换(6.48)就得到 $K_{ijk}$ 上的节点基函数.采用一维参考单元上Lagrange线性元节点基函数（参见表6.1）

$$
\Phi_ {0} (\xi) = 1 - \xi , \quad \Phi_ {1} (\xi) = \xi ,
$$

可得参考单元上的乘积型节点基函数 $\varPhi_l(\xi)\varPhi_m(\eta)\varPhi_n(\zeta),l,m,n=0,1.$ 从而得到 $K_{ijk}$ 上的节点基函数：

$$
\phi_ {i j k, l m n} (x, y) = \Phi_ {l} \left(\frac {x - x _ {i - 1}}{h _ {x i}}\right) \Phi_ {m} \left(\frac {y - y _ {j - 1}}{h _ {y j}}\right) \Phi_ {k} \left(\frac {z - z _ {k - 1}}{h _ {z k}}\right). \tag {6.49}
$$

显然，任一有限元函数 $v \in Q_{1}(K_{ijk})$ 都可以表示为

$$
v (x, y, z) = \sum_ {l, m, n = 0} ^ {1} v \left(x _ {i - 1 + l}, y _ {j - 1 + m}, z _ {k - 1 + n}\right) \phi_ {i j k, l m n} (x, y, z), (x, y, z) \in K _ {i j k}. \tag {6.50}
$$

易见，由此插值公式生成的试探函数整体属于 $C$ ，因此相应的三线性有限元空间 $U_{h}\subset H^{1}$

# 6.2.6 习题

1. 证明下列分段二次多项式基函数属于 $C^1([0,1])$ ，并作出其图形

$$
\varphi_ {i} ^ {(0)} (x) = \left\{ \begin{array}{l l} 0, & 0 \leqslant x \leqslant x _ {i - 1}, x _ {i + 1} \leqslant x \leqslant 1, \\ 2 \left(\frac {x - x _ {i - 1}}{h}\right) ^ {2}, & x _ {i - 1} \leqslant x \leqslant x _ {i} - 0. 5 h, \\ - 2 \left(\frac {x - x _ {i}}{h}\right) ^ {2} + 1, & x _ {i} - 0. 5 h \leqslant x \leqslant x _ {i} + 0. 5 h, \\ 2 \left(\frac {x - x _ {i + 1}}{h}\right) ^ {2}, & x _ {i} + 0. 5 h \leqslant x \leqslant x _ {i + 1}, \end{array} \right.
$$

$$
\varphi_ {i} ^ {(1)} (x) = \left\{ \begin{array}{l l} 0, & 0 \leqslant x \leqslant x _ {i - 1} \text {或} x _ {i + 1} \leqslant x \leqslant 1 \\ - \frac {1}{2 h} \left(x - x _ {i - 1}\right) ^ {2}, & x _ {i - 1} \leqslant x \leqslant x _ {i} - 0. 5 h, \\ x - x _ {i} + \frac {3}{2 h} \left(x - x _ {i}\right) ^ {2}, & x _ {i} - 0. 5 h \leqslant x \leqslant x _ {i}, \\ x - x _ {i} - \frac {3}{2 h} \left(x - x _ {i}\right) ^ {2}, & x _ {i} \leqslant x \leqslant x _ {i} + 0. 5 h, \\ \frac {1}{2 h} \left(x - x _ {i + 1}\right) ^ {2}, & x _ {i} + 0. 5 h \leqslant x \leqslant x _ {i + 1}. \end{array} \right.
$$

2. 试就 $f = 1$ 具体写出有限元方程 (6.24) (6.25).

3. 证明矩形网上双二次、双三次Lagrange有限元试探函数空间包含于 $H^1$

4. 证明积分公式 (6.40).

5. 设 $l$ 是 $xy$ 平面上的直线, 方程为 $ax + by + c = 0 (a^2 + b^2 = 1)$ , $\pmb{n}$ 是 $l$ 的单位法向量. 若多项式 $\pmb{p}$ 具性质

$$
\left. \frac {\mathrm {d} ^ {i} p}{\mathrm {d} \boldsymbol {n} ^ {i}} \right| _ {l} = 0, \quad 0 \leqslant i \leqslant k,
$$

则 $p_{\cdot}$ 可用 $(ax + by + c)^{k + 1}$ 整除，即有多项式 $q(x,y)$ 使

$$
p (x, y) = (a x + b y + c) ^ {k + 1} q (x, y).
$$

# 6.3 二阶椭圆型方程的有限元法

确定了网络剖分和单元形状函数后，试探函数空间 $U_{h}$ 也就定了.本节讨论二阶椭圆型方程的有限元法

设 $\Omega \subset \mathbb{R}^d$ 是有界多面体区域（二维时指多边形，一维时指线段），考虑如下的二阶椭圆型问题：

$$
- \nabla \cdot (\kappa (\boldsymbol {x}) \nabla u (\boldsymbol {x})) = f (\boldsymbol {x}), \quad \boldsymbol {x} \in \Omega , \tag {6.51}
$$

在 $\Gamma = \partial \Omega$ 上给出下列边值条件之一：

$$
u | _ {\Gamma} = 0 \quad (\text {第 一 边 值 条 件 或 D i r i c h l e t 边 值 条 件}), \tag {6.52}
$$

$$
\left. \kappa (\boldsymbol {x}) \nabla u (\boldsymbol {x}) \cdot \boldsymbol {n} \right| _ {\Gamma} = 0 \quad \text {(第 二 边 值 条 件 或 N e u m a n n 边 值 条 件)}, \tag {6.53}
$$

$$
\left. \left(\boldsymbol {\kappa} (\boldsymbol {x}) \nabla u (\boldsymbol {x}) \cdot \boldsymbol {n} + q u\right) \right| _ {\Gamma} = 0 \quad (\text {第 三 边 值 条 件 或 R o b i n 边 值 条 件}), \tag {6.54}
$$

其中 $\kappa (\pmb {x}) > 0,f(\pmb {x}),q(\pmb {x})$ 都是给定的连续函数.显然，当 $q = 0$ 时，第三边值条件就化为第二边值条件.边值条件也可以这样给：在 $\varGamma$ 的一部分满足一种边值条件，其余部分满足另一种边值条件.当然下面的讨论也可以容易地推广到非齐次边界条件的情形

# 6.3.1 有限元离散

这一小节, 我们基于椭圆型问题的变分形式给出其Lagrange线性有限元离散. 当然可以容易地推广到其他类型的有限元, 比如Lagrange高次元或Hermite元

设 $\mathcal{M}_h$ 为 $\Omega$ 的一个三角剖分. 引入Lagrange线性元空间：

$$
U _ {h} = \left\{v _ {h} \in H ^ {1} (\Omega): v _ {h} | _ {K} \in P _ {1} (K), \forall K \in \mathcal {M} _ {h} \right\}, \quad U _ {h} ^ {0} = U _ {h} \cap H _ {0} ^ {1} (\Omega).
$$

先考虑Dirichlet边值条件.由5.4.3小节的推导，利用第一Green公式(5.7)，易得问题（6.51）和（6.52）的变分形式为：求 $u\in H_0^1 (\varOmega)$ 使得

$$
a (u, v) = (f, v), \quad \forall v \in H _ {0} ^ {1} (\Omega), \tag {6.55}
$$

其中双线性形式

$$
a (u, v) = \int_ {\Omega} \kappa (\boldsymbol {x}) \nabla u \cdot \nabla v d \boldsymbol {x}. \tag {6.56}
$$

得问题 (6.51) 和 (6.52) 的线性有限元离散: 求 $u_h \in U_h^0$ 使得

$$
a \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in U _ {h} ^ {0}. \tag {6.57}
$$

再考虑Neumann边界条件.此时问题(6.51)和(6.53)的变分形式为：求 $u\in H^{1}(\varOmega)$

使得

$$
a (u, v) = (f, v), \quad \forall v \in H ^ {1} (\Omega), \tag {6.58}
$$

问题 (6.51) 和 (6.53) 的线性有限元离散为: 求 $u_{h} \in U_{h}$ 使得

$$
a \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in U _ {h}. \tag {6.59}
$$

需要说明的是, 注意到如果 $u_{h}$ 是 (6.59) 的解, 那么 $u_{h}$ 加上任意常数也是 (6.59) 的解, 即纯 Neumann 问题的有限元解不唯一. 我们可以对 $u_{h}$ 加适当条件保证唯一性, 比如要求 $\int_{\Omega} u_{h} \mathrm{d}\pmb{x} = 0$ . 另外, 取 $v_{h} = 1$ 可知为了保证解的存在性, $f$ 应该满足条件 $\int_{\Omega} f \mathrm{d}\pmb{x} = 0$ .

对 Robin 边值条件, 由

$$
- \int_ {\Omega} \nabla \cdot (\kappa \nabla u) v d x = \int_ {\Omega} \kappa \nabla u \cdot \nabla v d x - \int_ {\Gamma} \kappa \nabla u \cdot n v = \int_ {\Omega} \kappa \nabla u \cdot \nabla v d x + \int_ {\Gamma} q u v d s
$$

得问题 (6.51) 和 (6.54) 的变分形式为: 求 $u \in H^{1}(\Omega)$ 使得

$$
\tilde {a} (u, v) = (f, v), \quad \forall v \in H ^ {1} (\Omega), \tag {6.60}
$$

其中双线性形式

$$
\tilde {a} (u, v) = a (u, v) + \int_ {\Gamma} q u v \mathrm {d} s. \tag {6.61}
$$

问题 (6.51) 和 (6.54) 的线性有限元离散为: 求 $u_{h} \in U_{h}$ 使得

$$
\tilde {a} \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in U _ {h}. \tag {6.62}
$$

# 6.3.2 有限元方程组的形成

这一小节我们考虑如何形成有限元方程组的刚度矩阵和右端向量（也称为荷载向量).为了简单起见，以二维情形 $(d = 2)$ Lagrange线性有限元方法为例来讨论.

先考虑Neumann边值条件的情形.设 $\{P_j\}_{j = 1}^J$ 是网格 $\mathcal{M}_h$ 的节点的集合.设 $\{\phi_j\}_{j = 1}^J$ 线性元空间 $U_{h}$ 的节点基，满足 $\phi_i(P_j) = \delta_{ij},i,j = 1,2,\dots ,J.$ 令

$$
u _ {h} = u _ {1} \phi_ {1} + u _ {2} \phi_ {2} + \dots + u _ {J} \phi_ {J}, \quad \text {其 中} \quad u _ {j} = u _ {h} (P _ {j}),
$$

在 (6.59) 中取 $v_{h} = \phi_{i}$ , 得有限元方程组:

$$
a \left(\phi_ {1}, \phi_ {i}\right) u _ {1} + a \left(\phi_ {2}, \phi_ {i}\right) u _ {2} + \dots + a \left(\phi_ {J}, \phi_ {i}\right) u _ {J} = \left(f, \phi_ {i}\right), \quad i = 1, 2, \dots , J.
$$

简记 $a_{ij} = a(\phi_j,\phi_i),f_i = (f,\phi_i)$ ，可将有限元方程组写为矩阵形式：

$$
\begin{array}{r} {\pmb {A} \pmb {U} = \pmb {F}, \quad \text {其 中} \quad \pmb {A} = \left(a _ {i j}\right) _ {J \times J}, \quad \pmb {U} = \left(u _ {i}\right) _ {J \times 1}, \quad \pmb {F} = \left(f _ {i}\right) _ {J \times 1}.} \end{array} \tag {6.63}
$$

首先考虑刚度矩阵 $\mathbf{A}$ 的形成. 同前面一维模型问题, 我们先形成单元刚度矩阵, 再

装配总刚度矩阵. 我们有

$$
a _ {i j} = \sum_ {K \in \mathcal {M} _ {h}} \int_ {K} \kappa (\boldsymbol {x}) \nabla \phi_ {j} \cdot \nabla \phi_ {i} d \boldsymbol {x}. \tag {6.64}
$$

显然如果节点 $P_{i}$ 和 $P_{j}$ 不相邻, 则 $a_{ij} = 0$ , 即 $\mathbf{A}$ 是稀疏矩阵.

我们在每个单元 $K$ 上先计算形如 $\int_{K}\kappa (\pmb {x})\nabla \phi_{j}\cdot \nabla \phi_{i}\mathrm{d}\pmb{x}$ 的积分，再加起来得到 $a_{ij}$ 既然在每个单元 $K$ ，节点基函数的限制如果非零就是某个重心坐标函数 $L_{p},p = 1,2,3.$ 因此只需计算如下的 $3\times 3$ 矩阵

$$
\boldsymbol {A} ^ {K}: \quad a _ {p q} ^ {K} = \int_ {K} \kappa (\boldsymbol {x}) \nabla L _ {q} \cdot \nabla L _ {p} \mathrm {d} \boldsymbol {x}, \quad p, q = 1, 2, 3. \tag {6.65}
$$

这里 $\pmb{A}^{K}$ 称为单元刚度矩阵

为了找出单元刚度矩阵和总刚度矩阵之间的关系. 定义 $K_{p}$ $(p = 1,2,3)$ 为单元 $K$ 的第 $p$ 个顶点的总体编号，则 $\phi_{K_p}\big|_K = L_p$ ，且总刚度矩阵可以由单元刚度矩阵装配起来：

$$
\sum_ {\substack {K, p, q \\ K _ {p} = i, K _ {q} = j}} a _ {p q} ^ {K} = a _ {i j}. \tag{6.66}
$$

然后考虑右端荷载向量 $\pmb{F}$ 的形成. 我们有

$$
f _ {i} = \sum_ {K \in \mathcal {M} _ {h}} \int_ {K} f \phi_ {i} \mathrm {d} x. \tag {6.67}
$$

同样先计算如下单元荷载向量：

$$
\boldsymbol {F} ^ {K}: \quad f _ {p} ^ {K} = \int_ {K} f L _ {p} \mathrm {d} \boldsymbol {x}, \quad p = 1, 2, 3. \tag {6.68}
$$

再按以下公式装配总荷载向量 $\pmb{F}$

$$
\sum_ {K, p, K _ {p} = i} f _ {p} ^ {K} = f _ {i}. \tag {6.69}
$$

例如，考虑如图6.13的三角网中的两个相邻节点 $P_{i}$ 和 $P_{j}$ . 单元 $K^{\mathrm{I}}, K^{\mathrm{II}}, \dots, K^{\mathrm{VI}}$ 的顶点的局部标号如图6.13所示，则

$$
K _ {2} ^ {\mathrm {I}} = K _ {1} ^ {\mathrm {I I}} = K _ {3} ^ {\mathrm {I I I}} = K _ {2} ^ {\mathrm {I V}} = K _ {1} ^ {\mathrm {V}} = K _ {3} ^ {\mathrm {V I}} = i, \quad K _ {3} ^ {\mathrm {I}} = K _ {2} ^ {\mathrm {V I}} = j.
$$

由(6.66)和(6.69)得

$$
\begin{array}{l} a _ {i j} = a _ {2 3} ^ {K ^ {\mathrm {I}}} + a _ {3 2} ^ {K ^ {\mathrm {V I}}}, \\ a _ {i i} = a _ {2 2} ^ {K ^ {\mathrm {I}}} + a _ {1 1} ^ {K ^ {\mathrm {I I}}} + a _ {3 3} ^ {K ^ {\mathrm {I I I}}} + a _ {2 2} ^ {K ^ {\mathrm {I V}}} + a _ {1 1} ^ {K ^ {\mathrm {V}}} + a _ {3 3} ^ {K ^ {\mathrm {V I}}}, \\ f _ {i} = f _ {2} ^ {K ^ {\mathrm {I}}} + f _ {1} ^ {K ^ {\mathrm {I I}}} + f _ {3} ^ {K ^ {\mathrm {I I I}}} + f _ {2} ^ {K ^ {\mathrm {I V}}} + f _ {1} ^ {K ^ {\mathrm {V}}} + f _ {3} ^ {K ^ {\mathrm {V I}}}. \\ \end{array}
$$

当然，实际编程装配总刚度矩阵时，并不是对其中元素的下标循环，即不是对每组 $i, j$ 找哪些单元刚度矩阵的元素对 $a_{ij}$ 有贡献，而是对单元和单元刚度矩阵元素的下标循环，

![](images/0ce2ae0d8dd9f78cadee8182753a66823cbacd5f9efc4fa45bc81e4ab0b3cc17.jpg)  
图6.13 整体与局部编号

利用公式 (6.66) 和 (6.69), 把相应的元素累加到总刚度矩阵对应的位置. 由单刚合成总刚 $A$ 的伪代码见算法 6.1.

# 算法6.1 装配总刚矩阵

令 $a_{ij} = 0, i,j = 1,2,\dots ,J.$

for $K\in \mathcal{M}_h$ $p,q = 1,2,3$ do

$$
a _ {K _ {p} K _ {q}} \leftarrow a _ {K _ {p} K _ {q}} + a _ {p q} ^ {K}
$$

end for

同样，可按算法6.2装配右端向量 $\pmb{F}$

# 算法6.2 装配右端向量

令 $f_{i} = 0, i = 1,2,\dots ,J.$

for $K\in \mathcal{M}_h$ $p = 1,2,3$ do

$$
f _ {K _ {p}} \leftarrow f _ {K _ {p}} + f _ {p} ^ {K}
$$

end for

下面考虑Dirichlet边值条件情形.不妨设 $\{P_j\}_{j = 1}^{J_0}$ 是所有内部节点的集合，即 $P_{j}\in$ $\partial \Omega ,j = J_0 + 1,J_0 + 2,\dots ,J$ 是边界节点.显然 $\{\phi_j\}_{j = 1}^{J_0}$ 组成了 $U_h^0$ 的节点基.记 $u_{h} =$ （20 $u_{1}^{0}\phi_{1} + u_{2}^{0}\phi_{2} + \dots +u_{J_{0}}^{0}\phi_{J_{0}}$ 可得有限元方法(6.57）的方程组：

$\pmb{A}^{0}\pmb{U}^{0} = \pmb{F}^{0}$ ，其中 $\pmb {A}^0 = \left(a_{ij}\right)_{J_0\times J_0},\quad \pmb {U}^0 = \left(u_i^0\right)_{J_0\times 1},\quad \pmb {F}^0 = \left(f_i\right)_{J_0\times 1}.$ (6.70)

显然按算法6.1和6.2得到 $\pmb{A}$ 和 $\pmb{F}$ 之后，消掉 $\pmb{A}$ 的第 $J_0 + 1, J_0 + 2, \dots, J$ 行和列及 $\pmb{F}$ 的这些行之后，就可以得到Dirichlet边值问题(6.51)(6.52)的线性有限元的刚度矩阵 $A^0$ 和右端 $F^0$ 了.

对于 Robin 边值条件情形, 易知有限元方法 (6.57) 的方程组可表示为

$(A + Q)U = F,$ 其中 $\mathbf{A},\mathbf{F}$ 见(6.63)， $Q = \left(q_{ij}\right)_{J\times J},q_{ij} = \int_{\Gamma}q\phi_{j}\phi_{i}\mathrm{d}s.$ (6.71)

显然, 只需再装配矩阵 $Q$ 即可. 为此, 记 $\mathcal{E}_h^B$ 为边界 $\Gamma$ 上单元的边的集合. 对任一 $e \in \mathcal{E}_h^B$ , 设其第 $p$ 个顶点对应的重心坐标函数为 $L_p$ , 其总体编号记为 $e_p, p = 1,2$ . 显然 $\phi_{e_p}|_e = L_p$ . 先计算每个 $e$ 上的矩阵

$$
Q ^ {e}: \quad q _ {p q} ^ {e} = \int_ {e} q L _ {q} L _ {p} \mathrm {d} s, \quad p, q = 1, 2, \tag {6.72}
$$

就可以类似于前面的算法装配边值条件矩阵 $Q$

下面考虑积分 (6.65) 和 (6.68) 的计算. 对于复杂系数, 只能采用近似计算. 记 $S$ 为 $\triangle$ 的面积. 为了便于读者查阅, 在表 6.2 中列出了三角形 $\triangle(1,2,3)$ 上的 Gauss 求积公式

$$
\int_ {\triangle} f \mathrm {d} x \sim S \sum_ {i} W _ {i} f \left(L _ {1} ^ {(i)}, L _ {2} ^ {(i)}, L _ {3} ^ {(i)}\right)
$$

的求积节点的面积坐标 $\left(L_1^{(i)}, L_2^{(i)}, L_3^{(i)}\right)$ 和求积系数 $W_i$ . 公式关于面积坐标 $(L_1, L_2, L_3)$ 是对称分布的. 若表中出现的 $L_1, L_2, L_3$ 互异，则经过置换应有六个求积节点. 若表中出现的 $L_1, L_2, L_3$ 有两个互异，则经置换后应有三个求积节点. 点 $\left(\frac{1}{3}, \frac{1}{3}, \frac{1}{3}\right)$ 是重心，若取作求积节点，则只出现一次 [29].

表 6.2 求积公式节点及系数表  

<table><tr><td>\( {W}_{i} \)</td><td>\( {L}_{1} \)</td><td>\( {L}_{2} \)</td><td>\( {L}_{3} \)</td><td>重数</td></tr><tr><td></td><td>3 点公式</td><td>2 阶精度</td><td></td><td></td></tr><tr><td>0.333 333 333 333 333</td><td>0.666 666 666 666 667</td><td>0.166 666 666 666 667</td><td>0.166 666 666 666 667</td><td>3</td></tr><tr><td></td><td>3 点公式</td><td>2 阶精度</td><td></td><td></td></tr><tr><td>0.333 333 333 333 333</td><td>0.500 000 000 000 000</td><td>0.500 000 000 000 000</td><td>0.000 000 000 000 000</td><td>3</td></tr><tr><td></td><td>4 点公式</td><td>3 阶精度</td><td></td><td></td></tr><tr><td>-0.562 500 000 000 000</td><td>0.333 333 333 333 333</td><td>0.333 333 333 333 333</td><td>0.333 333 333 333 333</td><td>1</td></tr><tr><td>0.520 833 333 333 333</td><td>0.600 000 000 000 000</td><td>0.200 000 000 000 000</td><td>0.200 000 000 000 000</td><td>3</td></tr><tr><td></td><td>6 点公式</td><td>3 阶精度</td><td></td><td></td></tr><tr><td>0.166 666 666 666 667</td><td>0.659 027 622 374 092</td><td>0.231 933 368 553 031</td><td>0.109 039 009 072 877</td><td>6</td></tr><tr><td></td><td>6 点公式</td><td>4 阶精度</td><td></td><td></td></tr><tr><td>0.109 951 743 655 322</td><td>0.816 847 572 980 459</td><td>0.091 576 213 509 771</td><td>0.091 576 213 509 771</td><td>3</td></tr><tr><td>0.223 381 589 678 011</td><td>0.108 103 018 168 070</td><td>0.445 948 490 915 965</td><td>0.445 948 490 915 965</td><td>3</td></tr><tr><td></td><td>7 点公式</td><td>4 阶精度</td><td></td><td></td></tr><tr><td>0.375 000 000 000 000</td><td>0.333 333 333 333 333</td><td>0.333 333 333 333 333</td><td>0.333 333 333 333 333</td><td>1</td></tr><tr><td>0.104 166 666 666 667</td><td>0.736 712 498 968 435</td><td>0.237 932 366 472 434</td><td>0.025 355 134 551 932</td><td>6</td></tr><tr><td></td><td>7 点公式</td><td>5 阶精度</td><td></td><td></td></tr><tr><td>0.225 030 000 330 000</td><td>0.333 333 333 333 333</td><td>0.333 333 333 333 333</td><td>0.333 333 333 333 333</td><td>1</td></tr><tr><td>0.125 939 180 544 827</td><td>0.797 426 985 353 087</td><td>0.101 286 507 323 456</td><td>0.101 286 507 323 456</td><td>3</td></tr><tr><td>0.132 394 152 788 506</td><td>0.470 142 064 105 115</td><td>0.470 142 064 105 115</td><td>0.059 715 871 789 770</td><td>3</td></tr><tr><td></td><td>9 点公式</td><td>5 阶精度</td><td></td><td></td></tr><tr><td>0.205 950 504 760 887</td><td>0.124 649 503 233 232</td><td>0.437 525 248 383 384</td><td>0.437 525 248 383 384</td><td>3</td></tr><tr><td>0.063 691 414 286 223</td><td>0.797 112 651 860 071</td><td>0.165 409 927 389 841</td><td>0.037 477 420 750 088</td><td>6</td></tr><tr><td></td><td>12 点公式</td><td>6 阶精度</td><td></td><td></td></tr><tr><td>0.050 844 906 370 207</td><td>0.873 821 971 016 996</td><td>0.063 089 011 491 502</td><td>0.063 089 011 491 502</td><td>3</td></tr><tr><td>0.116 786 275 726 379</td><td>0.501 426 509 658 179</td><td>0.249 286 745 170 910</td><td>0.249 286 745 170 910</td><td>3</td></tr><tr><td>0.082 851 075 618 374</td><td>0.636 502 499 121 399</td><td>0.310 352 451 033 785</td><td>0.053 145 049 844 816</td><td>6</td></tr><tr><td></td><td>13 点公式</td><td>7 阶精度</td><td></td><td></td></tr><tr><td>-0.149 570 044 467 670</td><td>0.333 333 333 333 333</td><td>0.333 333 333 333 333</td><td>0.333 333 333 333 333</td><td>1</td></tr><tr><td>0.175 615 257 433 204</td><td>0.479 308 067 841 923</td><td>0.260 345 966 079 038</td><td>0.260 345 966 079 038</td><td>3</td></tr><tr><td>0.053 347 235 608 839</td><td>0.869 739 794 195 568</td><td>0.065 130 102 902 216</td><td>0.065 130 102 902 216</td><td>3</td></tr><tr><td>0.077 113 760 890 257</td><td>0.638 444 188 569 809</td><td>0.312 865 496 004 875</td><td>0.048 690 315 425 316</td><td>6</td></tr></table>

# 6.3.3 习题

1. 给出由 (6.72) 中的 $Q^{e}$ 装配 (6.71) 中的矩阵 $Q$ 的伪代码  
2. 对三维Dirichlet边值条件情形, 给出装配线性有限元法刚度矩阵的算法  
3. (实习题) 用线性有限元法求下列问题的数值解：

$$
\left\{ \begin{array}{l l} \Delta u = - 2, & - 1 <   x, y <   1, \\ u (x, - 1) = u (x, 1) = 0, & - 1 <   x <   1, \\ u _ {x} (- 1, y) = 1, u _ {x} (1, y) = 0, & - 1 <   y <   1 \end{array} \right.
$$

(精确到小数点后第6位).

# *6.4 有限元法的收敛性理论

本节考虑求解二维、三维椭圆型方程的有限元法的收敛性.类似于6.1节的一维情形，我们将有限元解的误差估计转化为插值误差估计，并利用所谓尺度变换的技巧（scalingargument）来估计插值误差.

# 6.4.1 插值理论

本小节给出Lagrange型有限元的插值误差估计，Hermite型有限元的插值理论可类似推导，供感兴趣的同学自行推导.

一些辅助结果

设 $\Omega \subset \mathbb{R}^d$ . $P_k(\Omega)$ 是 $\Omega$ 上次数 $\leqslant k$ 的多项式集合. 首先我们引入Bramble-Hilbert引理, 将被用来在参考单元上估计插值误差.

引理6.2(Bramble-Hilbert) 设 $\Omega \subset \mathbb{R}^d$ 是有界多面体区域, $m \geqslant 1$ 是整数, $\Pi: H^{m}(\Omega) \mapsto H^{m}(\Omega)$ 是一个有界线性算子, 且满足 $\Pi p = p, \forall p \in P_{m-1}(\Omega)$ . 则存在常数 $C = C(m, \Omega)$ 使得

$$
\| v - \Pi v \| _ {H ^ {m} (\Omega)} \leqslant C | v | _ {H ^ {m} (\Omega)}, \quad \forall v \in H ^ {m} (\Omega).
$$

证明 与Poincaré不等式的证明类似(见定理5.4).记 $\varPi_{m}:H^{m}(\varOmega)\mapsto P_{m-1}(\varOmega)$ 为按 $H^{m}(\varOmega)$ 内积的正交投影算子，则

$$
\| v - \Pi v \| _ {H ^ {m} (\Omega)} = \| v - \Pi_ {m} v - \Pi (v - \Pi_ {m} v) \| _ {H ^ {m} (\Omega)} \leqslant C \| v - \Pi_ {m} v \| _ {H ^ {m} (\Omega)}.
$$

引入空间

$$
V = \{v \in H ^ {m} (\Omega): \Pi_ {m} (v) = 0 \}.
$$

则只需证明

$$
\left\| v \right\| _ {H ^ {m} (\Omega)} \leqslant C | v | _ {H ^ {m} (\Omega)}, \quad \forall v \in V.
$$

用反证法. 假设上面不等式不成立, 则存在序列 $\{v_{n}\} \subset V$ 使得

$$
\left\| v _ {n} \right\| _ {H ^ {m} (\Omega)} = 1, \quad \left| v _ {n} \right| _ {H ^ {m} (\Omega)} \leqslant \frac {1}{n}.
$$

由 $H^{m}(\varOmega)\hookrightarrow\hookrightarrow H^{m-1}(\varOmega)$ （见例5.3(ii))，知存在 $H^{m-1}(\varOmega)$ 中收敛的子序列，仍记为 $\{v_n\}$ ，则 $\{v_n\}$ 是空间 $H^{m-1}(\varOmega)$ 中的Cauchy序列。再由 $|v_{n}|_{H^{m}(\varOmega)}\to 0$ ，知 $\{v_{n}\}$ 是空间 $H^{m}(\varOmega)$ 中的Cauchy序列，所以 $\{v_{n}\}$ 在 $H^{m}(\varOmega)$ 中收敛，其极限记为 $v$ ，则

$$
\left| v \right| _ {H ^ {m} (\Omega)} = 0, \quad \Pi_ {m} v = 0.
$$

由引理5.3知 $v\in P_{m - 1}(\varOmega)$ ，从而 $v = \Pi_{m}v = 0$ ，但 $\| v\|_{H^m (\varOmega)} = 1$ ，矛盾

事实上，Poincaré-Friedrichs 不等式可以作为 Bramble-Hilbert 引理的直接推论，其证明留作课后习题.

我们称两个单元 $K$ 和 $\hat{K} \subset \mathbb{R}^d$ 是仿射等价，如果存在可逆的仿射变换

$$
F: \hat {K} \to K, \quad F \hat {x} = B \hat {x} + b.
$$

对任意一个 $K$ 上的函数 $v(x)$ ，记 $\hat{v} = v\circ F.$ 显然 $\hat{v} (\hat{x}) = v(x).$ 下面两个引理将被用于尺度变换技巧.

引理6.3 设 $K$ 和 $\hat{K} \subset \mathbb{R}^d$ 是仿射等价的, 则 $\forall v \in H^{m}(K)$ , 有 $\hat{v} \in H^{m}(\hat{K})$ , 且存在常数 $C = d^{\frac{m}{2}}$ 使得

$$
\left| \hat {v} \right| _ {H ^ {m} (\hat {K})} \leqslant C \| B \| ^ {m} | \det  B | ^ {- \frac {1}{2}} | v | _ {H ^ {m} (K)},
$$

$$
\left| v \right| _ {H ^ {m} (K)} \leqslant C \left\| B ^ {- 1} \right\| ^ {m} | \det  B | ^ {\frac {1}{2}} | \hat {v} | _ {H ^ {m} (\hat {K})}.
$$

这里 $\| \cdot \|$ 为 $\mathbb{R}^d$ 中Euclid范数的诱导范数

证明 只证第一个不等式, 另一个同理可证. 由 $H^{m}$ 半范数的定义,

$$
| \hat {v} | _ {H ^ {m} (\hat {K})} ^ {2} = \int_ {\hat {K}} \sum_ {1 \leqslant i _ {1}, i _ {2}, \dots , i _ {m} \leqslant d} \left| \hat {\partial} _ {i _ {1}} \hat {\partial} _ {i _ {2}} \dots \hat {\partial} _ {i _ {m}} \hat {v} (\hat {x}) \right| ^ {2} \mathrm {d} \hat {x},
$$

$$
| v | _ {H ^ {m} (K)} ^ {2} = \int_ {K} \sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} \left| \partial_ {j _ {1}} \partial_ {j _ {2}} \dots \partial_ {j _ {m}} v (x) \right| ^ {2} \mathrm {d} x.
$$

由链式法则可得 $\hat{\partial}_i\hat{v} (\hat{x}) = \frac{\partial}{\partial\hat{x}_i}\hat{v} (\hat{x}) = \sum_{j = 1}^d\frac{\partial}{\partial x_j} v(x)\frac{\partial x_j}{\partial\hat{x}_i} = \sum_{j = 1}^d b_{ji}\partial_jv(x),$ 从而

$$
\hat {\partial} _ {i _ {1}} \hat {\partial} _ {i _ {2}} \dots \hat {\partial} _ {i _ {m}} \hat {v} (\hat {x}) = \sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} b _ {j _ {1} i _ {1}} b _ {j _ {2} i _ {2}} \dots b _ {j _ {m} i _ {m}} \partial_ {j _ {1}} \partial_ {j _ {2}} \dots \partial_ {j _ {m}} v (x).
$$

记 $d$ 阶单位矩阵的第 $i$ 列为 $e_i$ . 由Cauchy-Schwarz不等式得

$$
\begin{array}{l} \left| \hat {\partial} _ {i _ {1}} \hat {\partial} _ {i _ {2}} \dots \hat {\partial} _ {i _ {m}} \hat {v} (\hat {x}) \right| ^ {2} \\ \leqslant \left(\sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} | b _ {j _ {1} i _ {1}} b _ {j _ {2} i _ {2}} \dots b _ {j _ {m} i _ {m}} | ^ {2}\right) \sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} | \partial_ {j _ {1}} \partial_ {j _ {2}} \dots \partial_ {j _ {m}} v (x) | ^ {2} \\ = \left(\sum_ {j _ {1} = 1} ^ {d} \left| b _ {j _ {1} i _ {1}} \right| ^ {2} \sum_ {j _ {2} = 1} ^ {d} \left| b _ {j _ {2} i _ {2}} \right| ^ {2} \dots \sum_ {j _ {m} = 1} ^ {d} \left| b _ {j _ {m} i _ {m}} \right| ^ {2}\right) \sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} \left| \partial_ {j _ {1}} \partial_ {j _ {2}} \dots \partial_ {j _ {m}} v (x) \right| ^ {2} \\ = \left| B e _ {i _ {1}} \right| ^ {2} \left| B e _ {i _ {2}} \right| ^ {2} \dots \left| B e _ {i _ {m}} \right| ^ {2} \sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} \left| \partial_ {j _ {1}} \partial_ {j _ {2}} \dots \partial_ {j _ {m}} v (x) \right| ^ {2} \\ \leqslant \| B \| ^ {2 m} \sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} | \partial_ {j _ {1}} \partial_ {j _ {2}} \dots \partial_ {j _ {m}} v (x) | ^ {2}. \\ \end{array}
$$

最后两边积分并利用积分变换公式得

$$
\begin{array}{l} \left| \hat {v} \right| _ {H ^ {m} (\hat {K})} ^ {2} \leqslant d ^ {m} \| B \| ^ {2 m} \int_ {\hat {K}} \sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} \left| \partial_ {j _ {1}} \partial_ {j _ {2}} \dots \partial_ {j _ {m}} v (x) \right| ^ {2} \mathrm {d} \hat {x} \\ = d ^ {m} \| B \| ^ {2 m} \int_ {K} \sum_ {1 \leqslant j _ {1}, j _ {2}, \dots , j _ {m} \leqslant d} | \partial_ {j _ {1}} \partial_ {j _ {2}} \dots \partial_ {j _ {m}} v (x) | ^ {2} | \det  B ^ {- 1} | \mathrm {d} x \\ = d ^ {m} \| B \| ^ {2 m} \left| \det  B ^ {- 1} \right| | v | _ {H ^ {m} (K)} ^ {2}. \\ \end{array}
$$

两边开方得证

引理6.4 设 $K$ 和 $\hat{K}$ 是仿射等价的。记 $h_K = \mathrm{diam}(K), h_{\hat{K}} = \mathrm{diam}(\hat{K}), \rho_K$ 和 $\rho_{\hat{K}}$ 分别是 $K$ 和 $\hat{K}$ 的内接球直径，则如下估计成立

$$
\| B \| \leqslant \frac {h _ {K}}{\rho_ {\hat {K}}}, \quad \| B ^ {- 1} \| \leqslant \frac {h _ {\hat {K}}}{\rho_ {K}}, \quad | \det  B | = \frac {| K |}{| \hat {K} |}. \tag {6.73}
$$

证明 显然

$$
\| B\| = \frac{1}{\rho_{\hat{K}}}\sup_{|\xi | = \rho_{\hat{K}}}  |B\xi |.
$$

给定 $\xi \in \mathbb{R}^d$ 满足 $|\xi| = \rho_{\hat{K}}$ ，存在 $\hat{y}, \hat{z} \in \hat{K}$ 使得 $\hat{y} - \hat{z} = \xi$ （如图6.14），则 $F(\hat{y}), F(\hat{z}) \in K$ 且 $B\xi = F(\hat{y}) - F(\hat{z})$ 。从而 $|B\xi| \leqslant h_K$ 。（4.5）中第一个不等式得证。第二个不等式可由第一个得到。最后一个等式由积分变换公式（取被积函数为1）可得。

![](images/0bbaa4163d0a332adcf228638da37633a9f20870d71185540fcb4a99c6e7fca6.jpg)

![](images/ea0c920dd5b3994f1768750a5e7c014d791e98d389798a6dbf0b611ee8819035.jpg)  
图6.14 $K$ 和 $\hat{K}$ 之间的仿射变换

局部有限元插值

回忆6.2节，知一个单元上的一般的有限元插值的定义6.3具体到Lagrange型有限元，可改写为

定义6.5（局部Lagrange型有限元插值）给定Lagrange型有限元 $(K,\mathcal{P},\mathcal{N})$ ，设相应的插值节点为 $A_{1},A_{2},\dots ,A_{n}$ ，其中 $n = \dim \mathcal{P}$ ，则一个连续函数 $\pmb{\upsilon}$ 在单元 $K$ 上的局部有限元插值定义为

$$
I _ {K} v \in \mathcal {P}: \quad (I _ {K} v) (A _ {i}) = v (A _ {i}), \quad i = 1, 2, \dots , n.
$$

显然当 $v \in \mathcal{P}$ 时, $I_K v = v$ . 记 $\phi_i \in \mathcal{P}, i = 1,2,\dots,n$ 为节点基函数, 则 $I_K v = \sum_{i=1}^{n} v(A_i) \phi_i$ .

例6.3（Lagrange线性元插值）设 $K\subset \mathbb{R}^d$ 为单纯形 $A_{1}A_{2}\dots A_{d + 1}$ ，其上线性有限元空间 $P_{1}(K)$ 的节点基为 $\{L_i,i = 1,2,\dots ,d + 1\}$ ，则Lagrange插值为

$$
(I _ {K} v) (x) := \sum_ {i = 1} ^ {d + 1} v (A _ {i}) L _ {i} (x), \quad \forall v \in C (K).
$$

进行坐标变换, 将 $x$ 变为 $\hat{x} = (L_1(x), L_2(x), \dots, L_d(x))$ . 显然, $\hat{x}$ 由重心坐标的前 $d$ 个分量组成. 单纯形 $A_1A_2 \dots A_{d+1}$ 变为参考单元 $\hat{K}$ 即单纯形 $\hat{A}_1\hat{A}_2 \dots \hat{A}_{d+1}$ , 其中 $\hat{A}_1 = (1,0,\dots,0), \hat{A}_d = (0,0,\dots,1), \hat{A}_{d+1} = (0,0,\dots,0)$ . 函数 $v(x)$ 变为 $\hat{v}(\hat{x}) = v(x)$ , $v$ 在单元 $K$ 上的插值 $I_K v$ 变为 $\hat{v}(\hat{x})$ 在参考单元 $\hat{K}$ 上的插值 $I_{\hat{K}} \hat{v} \in P_1(\hat{K})$ , 满足:

$$
(I _ {\hat {K}} \hat {v}) (\hat {A} _ {i}) = \hat {v} (\hat {A} _ {i}), \quad i = 1, 2, \dots , d + 1.
$$

先考虑参考单元上的有限元插值估计

定理6.5 给定Lagrange型有限元 $(\hat{K},\hat{\mathcal{P}},\hat{\mathcal{N}})$ ，设

(i) $P_{m - 1}(\hat{K})\subset \hat{\mathcal{P}}\subset H^{m}(\hat{K})$   
(ii) $m - \frac{d}{2} >0$

则存在常数 $C = C(d,m,\hat{K})$ ，对 $0\leqslant i\leqslant m$ 有

$$
\left| \hat {v} - \hat {I} \hat {v} \right| _ {H ^ {i} (\hat {K})} \leqslant C \left| \hat {v} \right| _ {H ^ {m} (\hat {K})}, \quad \forall \hat {v} \in H ^ {m} (\hat {K}),
$$

其中 $\hat{I} = I_{\hat{K}}$ 是相应的局部有限元插值算子.

证明 首先, 由 Sobolev 嵌入定理 5.3, $H^{m}(\hat{K}) \hookrightarrow C(\hat{K})$ , 所以 $\forall \hat{v} \in H^{m}(\hat{K}), \hat{I}\hat{v}$ 有定义, 且 $\| \hat{I}\hat{v} \|_{H^{m}(\hat{K})} \leqslant C\|\hat{u}\|_{C(\hat{K})} \leqslant C\|\hat{u}\|_{H^{m}(\hat{K})}$ , 即 $\hat{I}$ 是由 $H^{m}(\hat{K})$ 到 $H^{m}(\hat{K})$ 的有界线性算子. 注意到 $I_{K}v = v, \forall v \in P_{m-1}(\hat{K})$ , 由 Bramble-Hilbert 引理得

$$
\| \hat {v} - \hat {I} \hat {v} \| _ {H ^ {m} (\hat {K})} \leqslant C | \hat {v} | _ {H ^ {m} (\hat {K})}, \quad \forall \hat {v} \in H ^ {m} (\hat {K}).
$$

证毕.

定义6.6（仿射插值等价）给定有限元 $(K,\mathcal{P},\mathcal{N})$ 和参考有限元 $(\hat{K},\hat{\mathcal{P}},\hat{\mathcal{N}})$ .相应的插值算子分别记为 $I = I_{K}$ 和 $\hat{I} = I_{\hat{K}}$ .如果存在可逆的仿射变换 $x = F(\hat{x}) = B\hat{x} +b$ 使得

(i) $K = F(\hat{K})$   
(ii) $\mathcal{P} = \{p:\hat{p}\in \hat{\mathcal{P}}\}$   
(iii) $\widehat{Iv} = \hat{I}\hat{v},$

那么称这两个有限元是仿射插值等价的.

显然例6.3中的Lagrange线性元满足仿射插值等价性质

定理6.6 设参考有限元 $(\hat{K},\hat{\mathcal{P}},\hat{\mathcal{N}})$ 满足定理6.5的条件，有限元 $(K,\mathcal{P},\mathcal{N})$ 仿射插值等价于 $(\hat{K},\hat{\mathcal{P}},\hat{\mathcal{N}})$ ，则存在 $C = C\left(d,m,\hat{K},\frac{h_K}{\rho_K}\right)$ ，对 $0\leqslant i\leqslant m$ 有

$$
\left| v - I _ {K} v \right| _ {H ^ {i} (K)} \leqslant C h _ {K} ^ {m - i} | v | _ {H ^ {m} (K)}, \quad \forall v \in H ^ {m} (K).
$$

证明 简记 $I = I_{K}$ ．由引理6.3—6.4及定理6.5得

$$
\left| v - I v \right| _ {H ^ {i} (K)} \leqslant C \left\| B ^ {- 1} \right\| ^ {i} | \det  B | ^ {\frac {1}{2}} \left| \widehat {v - I v} \right| _ {H ^ {i} (\hat {K})}
$$

$$
\begin{array}{l} = C \left\| B ^ {- 1} \right\| ^ {i} | \det  B | ^ {\frac {1}{2}} | \hat {v} - \hat {I} \hat {v} | _ {H ^ {i} (\hat {K})} \\ \leqslant C \left\| B ^ {- 1} \right\| ^ {i} | \det  B | ^ {\frac {1}{2}} | \hat {v} | _ {H ^ {m} (\hat {K})} \\ \leqslant C \left\| B ^ {- 1} \right\| ^ {i} \| B \| ^ {m} | v | _ {H ^ {m} (K)} \\ \leqslant C \left(\frac {h _ {\hat {K}}}{\rho_ {K}}\right) ^ {i} \left(\frac {h _ {K}}{\rho_ {\hat {K}}}\right) ^ {m} | v | _ {H ^ {m} (K)} \leqslant C \frac {h _ {\hat {K}} ^ {i}}{\rho_ {\hat {K}} ^ {m}} \left(\frac {h _ {K}}{\rho_ {K}}\right) ^ {i} h _ {K} ^ {m - i} | v | _ {H ^ {m} (K)}. \\ \end{array}
$$

得证.

例6.4（Lagrange线性元——局部插值误差）Lagrange线性元应用定理6.6只需验证定理6.5的条件.当 $d = 1$ 时， $m$ 可以取1或2;当 $d = 2,3$ 时， $m = 2.$ 故由定理6.6可得如下 $L^2$ 和 $H^{1}$ 插值误差估计

$$
\| v - I _ {K} v \| _ {L ^ {2} (K)} + h _ {K} \| v - I _ {K} v \| _ {H ^ {1} (K)} \leqslant \left\{ \begin{array}{l l} C h _ {K} | v | _ {H ^ {1} (K)}, & d = 1, \\ C h _ {K} ^ {2} | v | _ {H ^ {2} (K)}, & d = 1, 2, 3. \end{array} \right. \tag {6.74}
$$

整体有限元插值

定义6.7（有限元插值）设 $\mathcal{M}_h$ 为区域 $\Omega$ 的一个网格剖分.对任意单元 $K\in \mathcal{M}_h$ 给定有限元 $(K,\mathcal{P}_K,\mathcal{N}_K)$ 和相应的局部插值算子 $I_{K}$ ，则（整体）有限元插值算子 $I_{h}$ 定义为

$$
\left. \left(I _ {h} v\right) \right| _ {K} = I _ {K} \left(v \mid_ {K}\right), \quad \forall K \in \mathcal {M} _ {h}.
$$

定义6.8 如果一族网格 $\{\mathcal{M}_h\}$ 满足，存在常数 $\mu > 0$ 使得

$$
\frac {h _ {K}}{\rho_ {K}} \leqslant \mu , \quad \forall K \in \mathcal {M} _ {h},
$$

则称这族网格是正则(regular)的.如果存在常数 $\nu >0$ 使得

$$
\frac {h}{h _ {K}} \leqslant \nu , \quad \forall K \in \mathcal {M} _ {h}, \quad h := \max  _ {K \in \mathcal {M} _ {h}} h _ {K},
$$

则称这族网格是拟均匀（quasi-uniform）的.

由定理6.6可得如下插值误差估计

定理6.7 假设 $\{\mathcal{M}_h\}$ 是多面体区域 $\Omega \subset \mathbb{R}^d$ 的一族正则剖分。设参考有限元 $(\hat{K},\hat{\mathcal{P}},\hat{\mathcal{N}})$ 满足定理6.5的条件。对每一个单元 $K\in \mathcal{M}_h$ ，假设有限元 $(K,\mathcal{P}_K,\mathcal{N}_K)$ 仿射插值等价于 $(\hat{K},\hat{\mathcal{P}},\hat{\mathcal{N}})$ 。记 $h = \max_{K\in \mathcal{M}_h}h_K$ ，则对 $0\leqslant i\leqslant m$ ，存在常数 $C = C(\hat{K},d,m,\mu) > 0$ 使得

$$
\left(\sum_ {K \in \mathcal {M} _ {h}} \| v - I _ {h} v \| _ {H ^ {i} (K)} ^ {2}\right) ^ {\frac {1}{2}} \leqslant C h ^ {m - i} | v | _ {H ^ {m} (\Omega)}, \quad \forall v \in H ^ {m} (\Omega).
$$

例6.5 (Lagrange线性元的插值误差) 显然，对正则三角剖分上的Lagrange线性

有限元，有如下 $L^2$ 和 $H^{1}$ 插值误差估计

$$
\| v - I _ {h} v \| _ {L ^ {2} (\Omega)} + h \| v - I _ {h} v \| _ {H ^ {1} (\Omega)} \leqslant \left\{ \begin{array}{l l} C h | v | _ {H ^ {1} (\Omega)}, & d = 1, \\ C h ^ {2} | v | _ {H ^ {2} (\Omega)}, & d = 1, 2, 3. \end{array} \right. \tag {6.75}
$$

有限元逆估计

我们知道，对一般的函数，其导数值不一定可以被函数值来控制，但对有限元函数可以.

定理6.8 设单元 $K$ 仿射等价于参考单元 $\hat{K}$ 对 $m\geqslant 0$ 有

$$
\| \nabla v \| _ {L ^ {2} (K)} \leqslant C h _ {K} ^ {- 1} \| v \| _ {L ^ {2} (K)}, \quad \forall v \in P _ {m} (K), \tag {6.76}
$$

其中 $C$ 仅依赖于 $d,m,\hat{K}$ 及 $\frac{h_K}{\rho_K}$

证明 由引理6.3及有限维空间任意两个范数等价，我们有

$$
\begin{array}{l} \| \nabla v \| _ {L ^ {2} (K)} = | v | _ {H ^ {1} (K)} \leqslant C \left\| B _ {K} ^ {- 1} \right\| \left| \det  B _ {K} \right| ^ {\frac {1}{2}} | \hat {v} | _ {H ^ {1} (\hat {K})} \\ \leqslant C \left\| B _ {K} ^ {- 1} \right\| | \det  B _ {K} | ^ {\frac {1}{2}} \| \hat {v} \| _ {H ^ {1} (\hat {K})} \\ \leqslant C \left\| B _ {K} ^ {- 1} \right\| \left| \det  B _ {K} \right| ^ {\frac {1}{2}} \| \hat {v} \| _ {L ^ {2} (\hat {K})} \\ \leqslant C \left\| B _ {K} ^ {- 1} \right\| \| v \| _ {L ^ {2} (K)} \\ \leqslant C \frac {h _ {\hat {K}}}{\rho_ {K}} \| v \| _ {L ^ {2} (K)} = C h _ {\hat {K}} \frac {h _ {K}}{\rho_ {K}} h _ {K} ^ {- 1} \| v \| _ {L ^ {2} (K)}. \\ \end{array}
$$

得证.

局部迹不等式

定理6.9 设 $K$ 是 $\mathbb{R}^d$ 中直径为 $h_K$ 的多面体单元，仿射等价于 $\hat{K}$ ，则下面迹不等式成立：

$$
\begin{array}{l} \| v \| _ {L ^ {2} (\partial K)} \leqslant C \left(h _ {K} ^ {- \frac {1}{2}} \| v \| _ {L ^ {2} (K)} + \| v \| _ {L ^ {2} (K)} ^ {\frac {1}{2}} \| \nabla v \| _ {L ^ {2} (K)} ^ {\frac {1}{2}}\right) \\ \leqslant C \left(h _ {K} ^ {- \frac {1}{2}} \| v \| _ {L ^ {2} (K)} + h _ {K} ^ {\frac {1}{2}} \| \nabla v \| _ {L ^ {2} (K)}\right), \quad \forall v \in H ^ {1} (K), \tag {6.77} \\ \end{array}
$$

其中常数 $C$ 依赖于单元 $K$ 的正则性和 $\hat{K}$

证明设 $\hat{K}$ 的任一面 $\hat{e}$ 在仿射变换下变为 $K$ 的面 $e$ . 由定理5.5及引理6.3和6.4得

$$
\begin{array}{l} \| v \| _ {L ^ {2} (\partial K)} ^ {2} = \sum_ {e \subset \partial K} \int_ {e} v ^ {2} \mathrm {d} s = \sum_ {\hat {e} \subset \partial \hat {K}} \int_ {\hat {e}} \hat {v} ^ {2} \frac {| e |}{| \hat {e} |} \mathrm {d} s \leqslant \max  _ {e \subset \partial K} \frac {| e |}{| \hat {e} |} \| \hat {v} \| _ {L ^ {2} (\partial \hat {K})} ^ {2} \\ \leqslant C \max  _ {e \subset \partial K} \frac {| e |}{| \hat {e} |} \| \hat {v} \| _ {L ^ {2} (\hat {K})} \| \hat {v} \| _ {H ^ {1} (\hat {K})} \\ \end{array}
$$

$$
\begin{array}{l} \leqslant C \max  _ {e \subset \partial K} \frac {| e |}{| \hat {e} |} \left(\| \hat {v} \| _ {L ^ {2} (\hat {K})} ^ {2} + \| \hat {v} \| _ {L ^ {2} (\hat {K})} | \hat {v} | _ {H ^ {1} (\hat {K})}\right) \\ \leqslant C \max  _ {e \subset \partial K} \frac {| e |}{| \hat {e} |} | \det  B | ^ {- 1} \left(\| v \| _ {L ^ {2} (K)} ^ {2} + \| B \| \| v \| _ {L ^ {2} (K)} | v | _ {H ^ {1} (K)}\right) \\ \leqslant C \max  _ {e \subset \partial K} \frac {| e |}{| \hat {e} |} \frac {| \hat {K} |}{| K |} \left(\| v \| _ {L ^ {2} (K)} ^ {2} + \frac {h _ {K}}{\rho_ {\hat {K}}} \| v \| _ {L ^ {2} (K)} | v | _ {H ^ {1} (K)}\right) \\ \leqslant C \frac {h _ {K} ^ {d - 1}}{h _ {K} ^ {d}} \left(\| v \| _ {L ^ {2} (K)} ^ {2} + h _ {K} \| v \| _ {L ^ {2} (K)} | v | _ {H ^ {1} (K)}\right). \\ \end{array}
$$

得证.

# 6.4.2 误差估计

$H^{1}$ 误差估计

设 $\Omega$ 是 $\mathbb{R}^d$ 中的多边形区域, $\{\mathcal{M}_h\}$ 是 $\Omega$ 的一族三角剖分. 设 $U_h \subset H_0^1(\Omega)$ 是 $\mathcal{M}_h$ 上的协调线性有限元空间. 设 $u \in H_0^1(\Omega)$ 是下面变分问题的弱解:

$$
a (u, v) = \langle f, v \rangle , \quad \forall v \in H _ {0} ^ {1} (\Omega). \tag {6.78}
$$

$u_{h}\in U_{h}$ 是对应的有限元解：

$$
a \left(u _ {h}, v _ {h}\right) = \langle f, v _ {h} \rangle , \quad \forall v _ {h} \in U _ {h}. \tag {6.79}
$$

假设双线性形式 $a:H_0^1 (\varOmega)\times H_0^1 (\varOmega)\to \mathbb{R}$ 有界的和 $H_0^1 (\varOmega)$ -强制的：

$$
\left| a (u, v) \right| \leqslant \beta \| u \| _ {H ^ {1} (\Omega)} \| v \| _ {H ^ {1} (\Omega)}, \quad a (v, v) \geqslant \alpha \| v \| _ {H ^ {1} (\Omega)} ^ {2}, \quad \forall u, v \in H _ {0} ^ {1} (\Omega). \tag {6.80}
$$

定理6.10 假设解 $u\in H_0^1 (\varOmega)\cap H^2 (\varOmega)$ ，则存在与 $h$ 无关的常数 $C$ 使得

$$
\left\| u - u _ {h} \right\| _ {H ^ {1} (\Omega)} \leqslant C h \left| u \right| _ {H ^ {2} (\Omega)}.
$$

证明 由引理5.5及有限元插值误差估计(6.75)，

$$
\| u - u _ {h} \| _ {H ^ {1} (\Omega)} \leqslant C \inf  _ {v _ {h} \in U _ {h}} \| u - v _ {h} \| _ {H ^ {1} (\Omega)} \leqslant C \| u - I _ {h} u \| _ {H ^ {1} (\Omega)} \leqslant C h | u | _ {H ^ {2} (\Omega)}.
$$

证毕.

下面定理说明：即使(6.78)的解不属于 $H^2 (\Omega)$ ，有限元解仍然收敛

定理6.11

$$
\lim  _ {h \to 0} \| u - u _ {h} \| _ {H ^ {1} (\Omega)} = 0.
$$

证明 只需证明

$$
\lim  _ {h \rightarrow 0} \inf  _ {v _ {h} \in U _ {h}} \| u - v _ {h} \| _ {H ^ {1} (\Omega)} = 0.
$$

对任意 $\varepsilon > 0$ 及 $u \in H_0^1(\Omega)$ , 存在函数 $u_\varepsilon \in C_0^\infty(\Omega)$ 使得

$$
\| u - u _ {\varepsilon} \| _ {H ^ {1} (\Omega)} <   \frac {\varepsilon}{2}.
$$

另外，由插值误差估计 (6.75)，

$$
\left\| u _ {\varepsilon} - I _ {h} u _ {\varepsilon} \right\| _ {H ^ {1} (\Omega)} \leqslant C h \left| u _ {\varepsilon} \right| _ {H ^ {2} (\Omega)}.
$$

故存在 $h_\varepsilon > 0$ , 使得当 $0 < h < h_\varepsilon$ 时,

$$
\left\| u _ {\varepsilon} - I _ {h} u _ {\varepsilon} \right\| _ {H ^ {1} (\Omega)} <   \frac {\varepsilon}{2}.
$$

因此

$$
\inf  _ {v _ {h} \in U _ {h}} \| u - v _ {h} \| _ {H ^ {1} (\Omega)} \leqslant \| u - I _ {h} u _ {\varepsilon} \| _ {H ^ {1} (\Omega)} <   \frac {\varepsilon}{2} + \frac {\varepsilon}{2} = \varepsilon .
$$

证毕.

![](images/38d6740695954b02f7b1f187e227a6ad84f3e840a9cbefef57a6134c740b6a8b.jpg)

$L^2$ 误差估计

假设 (6.78) 的伴随问题在以下意义下是正则的: 对任意 $g \in L^{2}(\Omega)$ , 伴随问题

$$
a (v, \varphi_ {g}) = (v, g), \quad \forall v \in H _ {0} ^ {1} (\Omega)
$$

存在唯一解 $\varphi_g\in H^2 (\varOmega)\cap H_0^1 (\varOmega)$ ；并存在常数 $C$ 使得 $\| \varphi_{g}\|_{H^{2}(\varOmega)}\leqslant C\| g\|_{L^{2}(\varOmega)}$

定理6.12 假设问题(6.78)的解 $u\in H^{2}(\Omega)$ 且其伴随问题是正则的，则存在与 $h$ 无关的常数 $C$ 使得

$$
\left\| u - u _ {h} \right\| _ {L ^ {2} (\Omega)} \leqslant C h ^ {2} | u | _ {H ^ {2} (\Omega)}.
$$

证明 令 $g = u - u_{h},\varphi_{g}$ 为(4.10)的解，则

$$
\begin{array}{l} (u - u _ {h}, g) = a (u - u _ {h}, \varphi_ {g}) = a (u - u _ {h}, \varphi_ {g} - I _ {h} \varphi_ {g}) \\ \leqslant \beta \| u - u _ {h} \| _ {H ^ {1} (\Omega)} \| \varphi_ {g} - I _ {h} \varphi_ {g} \| _ {H ^ {1} (\Omega)} \\ \leqslant C h ^ {2} | u | _ {H ^ {2} (\Omega)} | \varphi_ {g} | _ {H ^ {2} (\Omega)} \\ \leqslant C h ^ {2} \| g \| _ {L ^ {2} (\Omega)} | u | _ {H ^ {2} (\Omega)}. \\ \end{array}
$$

证毕.

![](images/494cb39c1c317b36f84629cb232f0fff50e183bb7cdecdc6bb3d4efd9efa0f6e.jpg)

定理6.12的证明所用的技巧称为对偶论证或Aubin-Nitsche技巧

# 6.4.3 习题

1. 利用Bramble-Hilbert引理证明Poincaré-Friedrichs不等式  
2. 设单元 $K$ 与 $\hat{K} \subset \mathbb{R}^d$ 仿射等价, 证明局部的 Poincaré 不等式:

$$
\| v - v _ {K} \| _ {L ^ {2} (K)} \leqslant C h _ {K} \| \nabla v \| _ {L ^ {2} (K)}, \quad \forall v \in H ^ {1} (K),
$$

其中 $v_{K} = \frac{1}{|K|}\int_{K}v\mathrm{d}\pmb {x},C = C(\hat{K}).$

3. 给出矩形网格上双线性元插值的 $L^2$ 和 $H^1$ 误差估计

4. 证明 Dirichlet 边值问题 (6.51) (6.52) 对应的双线性形式 (6.56) 满足连续性和强制性条件 (6.80).  
5. 设 $U_h \subset H_0^1(\Omega)$ 是 $\mathcal{M}_h$ 上的二次有限元空间. 给出变分问题 (6.78) 的有限元离散 (6.79) 的 $H^1$ 和 $L^2$ 误差估计.

# 6.5 初边值问题的有限元法

Galerkin法以及由此发展起来的有限元法也可用于解初边值问题（非驻定问题），包括抛物型方程和双曲型方程.此时将时间变量 $t$ 看成参数，用虚功原理将初边值问题化成变分形式，然后用Galerkin有限元法求解

# 6.5.1 热传导方程

考虑热传导方程的初边值问题：

$$
\frac {\partial u}{\partial t} = \Delta u + f, \quad x \in \Omega \subset \mathbb {R} ^ {2}, t > 0, \tag {6.81}
$$

$$
u (\boldsymbol {x}; 0) = \psi (\boldsymbol {x}), \quad \boldsymbol {x} \in \Omega , \tag {6.82}
$$

$$
u | _ {\partial \Omega} = 0, \quad t > 0, \tag {6.83}
$$

其中 $\Delta u = \frac{\partial^2u}{\partial x^2} +\frac{\partial^2u}{\partial y^2},f = f(\pmb {x}),\Omega$ 是具分段光滑边界 $\partial \varOmega$ 的平面有界域.设对固定的 $t > 0,$ 解 $u(\pmb {x};t)$ 关于 $\pmb{x}$ 属于 $C^2 (\bar{\Omega})$ .以 $v\in H_0^1 (\varOmega)$ 乘(6.81）两端并积分，得

$$
\int_ {\Omega} \left(\frac {\partial u}{\partial t} - \Delta u - f\right) v d x = 0. \tag {6.84}
$$

利用Green公式和边值条件(6.83)，得

$$
\int_ {\Omega} \left(\frac {\partial u}{\partial t} v + \nabla u \cdot \nabla v - f v\right) \mathrm {d} x = 0. \tag {6.85}
$$

引进双线性形式和 $L^2 (\Omega)$ 内积：

$$
a (u, v) = \int_ {\Omega} \nabla u \cdot \nabla v \mathrm {d} x, \quad (f, v) = \int_ {\Omega} f v \mathrm {d} x.
$$

则初边值问题 (6.81)一(6.83)的变分形式为：求 $u(\cdot ;t)\in H_0^1 (\Omega)$ （视 $t$ 为参数)，满足

$$
\left(\frac {\partial u}{\partial t}, v\right) + a (u, v) = (f, v), \quad \forall v \in H _ {0} ^ {1} (\Omega), \tag {6.86}
$$

$$
u (\boldsymbol {x}; 0) = \psi (\boldsymbol {x}). \tag {6.87}
$$

并称如此的 $u$ 为初边值问题 (6.81) 的广义解

在 $H_0^1$ 中取一 $n$ 维子空间 $U_{h}$ ，所谓Galerkin法就是求含参数 $t$ 的函数 $u_{h}(\cdot ;t)\in$ $U_{h}(t\geqslant 0)$ ，满足

$$
\left(\frac {\partial u _ {h}}{\partial t}, v _ {h}\right) + a \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in U _ {h}, \tag {6.88}
$$

$$
(u _ {h} (\cdot ; 0), v _ {h}) = (\psi , v _ {h}), \quad \forall v _ {h} \in U _ {h}. \tag {6.89}
$$

在 $U_{h}$ 中取定一组基底 $\phi_1,\phi_2,\dots ,\phi_n$ ，将 $u_{h}$ 表示为

$$
u _ {h} (\boldsymbol {x}; t) = \sum_ {i = 1} ^ {n} \mu_ {i} (t) \phi_ {i} (\boldsymbol {x}),
$$

代到 (6.88), 并取 $v_{h} = \phi_{j}$ , 就得到关于 $\mu_1(t), \mu_2(t), \dots, \mu_n(t)$ 的常微分方程组:

$$
\sum_ {i = 1} ^ {n} \left(\phi_ {i}, \phi_ {j}\right) \frac {\mathrm {d} \mu_ {i}}{\mathrm {d} t} + \sum_ {i = 1} ^ {n} a \left(\phi_ {i}, \phi_ {j}\right) \mu_ {i} = (f, \phi_ {j}), \quad j = 1, 2, \dots , n. \tag {6.90}
$$

初值条件可按 (6.89) 或其他方法给出. 若按 (6.89) 取初值 $u_{h}(\pmb{x};0) = \sum_{i=1}^{n}\mu_{i}(0)\phi_{i}(\pmb{x})$ , 则 $\mu_{i}(0)$ 由下列方程组确定:

$$
\sum_ {i = 1} ^ {n} \left(\phi_ {i}, \phi_ {j}\right) \mu_ {i} (0) = (\psi , \phi_ {j}), \quad j = 1, 2, \dots , n. \tag {6.91}
$$

至此我们得到了半离散化 Galerkin 方程 (6.88) (6.89) 或 (6.90) (6.91). 若进一步对时间 $t$ 离散化, 就得到全离散化 Galerkin 法. 例如用向前差商 (向前 Euler 格式),

$$
\frac {1}{\tau} \left(u _ {h} ^ {n + 1} - u _ {h} ^ {n}, v _ {h}\right) + a \left(u _ {h} ^ {n}, v _ {h}\right) = (f, v _ {h}), \quad n = 0, 1, \dots , \tag {6.92}
$$

用向后差商 (向后 Euler 格式),

$$
\frac {1}{\tau} \left(u _ {h} ^ {n + 1} - u _ {h} ^ {n}, v _ {h}\right) + a \left(u _ {h} ^ {n + 1}, v _ {h}\right) = (f, v _ {h}), \quad n = 0, 1, \dots , \tag {6.93}
$$

或用 Crank-Nicolson 格式 (改进的 Euler 折线法)

$$
\frac {1}{\tau} \left(u _ {h} ^ {n + 1} - u _ {h} ^ {n}, v _ {h}\right) + \frac {1}{2} a \left(u _ {h} ^ {n + 1} + u _ {h} ^ {n}, v _ {h}\right) = (f, v _ {h}), \quad n = 0, 1, \dots , \tag {6.94}
$$

其中上标 $n$ 表示在 $t = t_{n} = n\tau$ 的近似

特别地，若对域 $\Omega$ 作三角剖分，并取 $U_h \subset H_0^1$ 为分片多项式函数空间，则前述

Galerkin法就是有限元法.以前各节关于有限元空间的构造和算法都可用到抛物型方程

注意全离散格式是由常微分方程组 (6.88) 或 (6.90) 离散化得到的, 所以需要讨论它的稳定性 (参看第 1 章). 作为例子, 我们证明格式 (6.94) 关于初值稳定. 在 (6.94) 中取 $v_{h} = u_{h}^{n + 1} - u_{h}^{n}, f = 0$ , 得

$$
\frac {1}{\tau} \left\| u _ {h} ^ {n + 1} - u _ {h} ^ {n} \right\| ^ {2} + \frac {1}{2} a \left(u _ {h} ^ {n + 1} + u _ {h} ^ {n}, u _ {h} ^ {n + 1} - u _ {h} ^ {n}\right) = 0, \tag {6.95}
$$

其中 $\| \cdot \|$ 表示 $L^2$ 范数. 利用 $a(\cdot, \cdot)$ 的对称性得

$$
a \left(u _ {h} ^ {n + 1} + u _ {h} ^ {n}, u _ {h} ^ {n + 1} - u _ {h} ^ {n}\right) = a \left(u _ {h} ^ {n + 1}, u _ {h} ^ {n + 1}\right) - a \left(u _ {h} ^ {n}, u _ {h} ^ {n}\right) = \left| u _ {h} ^ {n + 1} \right| _ {1} ^ {2} - \left| u _ {h} ^ {n} \right| _ {1} ^ {2}.
$$

其中 $|\cdot |_{1}$ 是 $H^{1}$ 半模.于是由(6.95）得

$$
\left| u _ {h} ^ {n + 1} \right| _ {1} \leqslant \left| u _ {h} ^ {n} \right| _ {1} \leqslant \dots \leqslant \left| u _ {h} ^ {0} \right| _ {1},
$$

由Poincaré-Friedrichs不等式知稳定性得证

# 6.5.2 波动方程

考虑波动方程的初边值问题：

$$
\frac {\partial^ {2} u}{\partial t ^ {2}} = \Delta u, \quad \boldsymbol {x} \in \Omega \subset \mathbb {R} ^ {2}, \quad t > 0, \tag {6.96}
$$

$$
u (\boldsymbol {x}; 0) = g _ {0} (\boldsymbol {x}), \quad u _ {t} (\boldsymbol {x}; 0) = g _ {1} (\boldsymbol {x}), \tag {6.97}
$$

$$
u | _ {\partial \Omega} = 0, \quad t > 0, \tag {6.98}
$$

其中 $\Delta u = \frac{\partial^2u}{\partial x^2} +\frac{\partial^2u}{\partial y^2},\Omega$ 是具分段光滑边界 $\partial \Omega$ 的平面有界域.与前述方法类似可得到(6.96)一(6.98）的变分形式

$$
\left(\frac {\partial^ {2} u}{\partial t ^ {2}}, v\right) + a (u, v) = 0, \quad \forall v \in H _ {0} ^ {1} (\Omega), \tag {6.99}
$$

$$
u (\boldsymbol {x}; 0) = g _ {0} (\boldsymbol {x}), \quad u _ {t} (\boldsymbol {x}; 0) = g _ {1} (\boldsymbol {x}), \tag {6.100}
$$

并称如此的解 $u$ 为初边值问题 (6.96)一(6.98)的广义解

在 $H_0^1$ 中取一 $n$ 维子空间 $U_{h}$ ，所谓Galerkin法就是求含参数 $t$ 的函数 $u_{h}(\cdot ;t)\in$ $U_{h}(t\geqslant 0)$ ，满足

$$
\left(\frac {\partial^ {2} u _ {h}}{\partial t ^ {2}}, v _ {h}\right) + a \left(u _ {h}, v _ {h}\right) = 0, \quad \forall v _ {h} \in U _ {h}, \tag {6.101}
$$

$$
(u _ {h} (\cdot ; 0), v _ {h}) = (g _ {0}, v _ {h}), \left(\frac {\partial u _ {h}}{\partial t} (\cdot ; 0), v _ {h}\right) = (g _ {1}, v _ {h}), \quad \forall v _ {h} \in U _ {h}. \tag {6.102}
$$

在 $U_{h}$ 中取定一组基底 $\varphi_1,\varphi_2,\dots ,\varphi_n$ 将 $u_{h}$ 表为

$$
u _ {h} (\boldsymbol {x}, t) = \sum_ {i = 1} ^ {n} \mu_ {i} (t) \varphi_ {i} (\boldsymbol {x}),
$$

代到 (6.101) (6.102), 并取 $v_{h} = \varphi_{j}$ , 就得到关于 $\mu_1, \mu_2, \dots, \mu_n$ 的常微分方程组:

$$
\sum_ {i = 1} ^ {n} \left(\varphi_ {i}, \varphi_ {j}\right) \frac {\mathrm {d} ^ {2} \mu_ {i}}{\mathrm {d} t ^ {2}} + \sum_ {i = 1} ^ {n} a \left(\varphi_ {i}, \varphi_ {j}\right) \mu_ {i} = 0, \quad j = 1, 2, \dots , n. \tag {6.103}
$$

初值条件为

$$
\mu_ {j} (0) = \alpha_ {j}, \quad \frac {\partial \mu_ {j} (0)}{\partial t} = \beta_ {j}. \tag {6.104}
$$

其中 $\mu_{j}(0)$ 和 $\frac{\partial\mu_j(0)}{\partial t}$ 按（6.101）（6.102）确定：

$$
\begin{array}{l} \sum_ {i = 1} ^ {n} \left(\varphi_ {i}, \varphi_ {j}\right) \mu_ {i} (0) = \left(g _ {0}, \varphi_ {j}\right), \tag {6.105} \\ \sum_ {i = 1} ^ {n} \left(\varphi_ {i}, \varphi_ {j}\right) \frac {\partial \mu_ {i} (0)}{\partial t} = \left(g _ {1}, \varphi_ {j}\right), \quad j = 1, 2, \dots , n. \\ \end{array}
$$

这样我们就得到半离散化 Galerkin 方程 (6.103) — (6.105). 这是一个二阶常微分方程组.

为得到全离散化 Galerkin 方程, 需用差分法对时间 $t$ 进一步离散. 设时间步长为

$$
\tau > 0, \quad t _ {n} = n \tau (n = 0, 1, \dots , N; N \tau = T),
$$

$u_h^n = u_h(x,y;t_n)$ .引进差商符号

$$
\begin{array}{l} \partial_ {t \bar {t}} u _ {h} ^ {n} = \frac {u _ {h} ^ {n + 1} - 2 u _ {h} ^ {n} + u _ {h} ^ {n - 1}}{\tau^ {2}}, \\ u _ {h} ^ {n, \frac {1}{4}} = \frac {u _ {h} ^ {n + 1} + 2 u _ {h} ^ {n} + u _ {h} ^ {n - 1}}{4}, \\ \end{array}
$$

则得全离散 Galerkin 方程, 例如, 显格式

$$
\left(\partial_ {t \bar {t}} u _ {h} ^ {n}, v _ {h}\right) + a \left(u _ {h} ^ {n}, v _ {h}\right) = 0, \quad \forall v _ {h} \in U _ {h}, \tag {6.106}
$$

隐格式

$$
\left(\partial_ {t \bar {t}} u _ {h} ^ {n}, v _ {h}\right) + a \left(u _ {h} ^ {n, \frac {1}{4}}, v _ {h}\right) = 0, \quad \forall v _ {h} \in U _ {h}. \tag {6.107}
$$

特别地, 若取 $U_h$ 为有限元空间, 则得 Galerkin 有限元法. 可以证明显格式 (6.106) 条件稳定, 隐格式 (6.107) 恒稳定 (参看 [3]).

（20 $1 = \frac{1}{4} x - 1.0 = \frac{1}{4} x - 1 <   0 <   \frac{1}{4} x$   
（2018·南京）如图，点 $P$ 是边长为 $3\mathrm{cm}$ 的圆周运动中，点 $Q$ 在圆心 $O$ 上，点 $R$ 在圆外，点 $S$ 在圆外，点 $T$ 在圆内，则点 $P,Q,R$ 之间的最大值是（ ）.

# 第7章

# 有限体积元法

有限体积元法也称广义差分法、控制体积法、有限体积法等. 最早由吉林大学李荣华教授提出 [19].

# 7.1 三角形网格上有限体积元法

考虑二阶椭圆型方程的边值问题：

$$
\left\{ \begin{array}{l} L u \equiv - \nabla \cdot (\kappa \nabla u) = f (x, y), \quad (x, y) \in \Omega , \\ u | _ {\partial \Omega} = 0, \end{array} \right. \tag {7.1}
$$

其中 $\Omega \subset \mathbb{R}^2$ 是一个多边形区域, $f \in L^{2}(\Omega)$ , 系数 $\kappa = (\kappa_{ij})_{i,j=1}^{2}$ 是一对称矩阵, 其元素 $\kappa_{ij}(x,y)$ 满足椭圆性条件, 即存在常数 $\gamma > 0$ , 使得

$$
\sum_ {i, j = 1} ^ {2} \kappa_ {i j} (x, y) \xi_ {i} \xi_ {j} \geqslant \gamma \sum_ {i = 1} ^ {2} \xi_ {i} ^ {2}
$$

对任意实向量 $(\xi_1,\xi_2)\in \mathbb{R}^2$ 和 $(x,y)\in \bar{\Omega}$ 成立

![](images/d17b1b83a0b3ae6659afacf2f2686ef3f8464fd9edd4885802613b7982167216.jpg)

人物简介

# 7.1.1 试探函数空间和检验函数空间

将 $\bar{\Omega}$ 分割成有限个小三角形之和，使不同三角形无重叠的内部区域，任一三角形的顶点不属于其他三角形的边的内部，且边界的角点都是三角形的顶点。每个三角形称为单元，三角形的顶点称为节点。所有单元构成 $\bar{\Omega}$ 的一个三角剖分，记为 $\mathcal{T}_h$ ，称为原始剖分。令 $h$ 表示所有三角形的最大边长。

现构造与 $\mathcal{T}_h$ 相关的对偶剖分 $\mathcal{T}_h^*$ ，常用的方法有两种：

(1) 重心对偶剖分. 取单元 $\triangle P_{0}P_{i}P_{i+1}$ 的重心 $Q_{i}$ 以及边中点 $M_{i}$ 为对偶剖分的节点, 则得重心对偶剖分, 如图7.1阴影部分所示.  
(2) 外心对偶剖分. 假设 $\tau_{h}$ 的每个单元的内角均不大于 $90^{\circ}$ , 取三角单元 $\triangle P_{0}P_{i}P_{i+1}$ 的外心 $Q_{i}$ 为对偶剖分的节点, 则得外心对偶剖分, 如图7.2阴影部分所示. 此时 $\overline{Q_{i}Q_{i+1}}$ 是 $\overline{P_{0}P_{i+1}}$ 的中垂线.

用 $\bar{\Omega}_h$ 表示剖分 $\mathcal{T}_h$ 的节点集合， $\dot{\Omega}_h = \bar{\Omega}_h \backslash \partial \Omega$ 表示内节点集合。 $\Omega_h^*$ 表示对偶剖分 $\mathcal{T}_h^*$ 的节点 $Q$ 的集合。对 $Q \in \Omega_h^*$ ，用 $K_Q$ 表示含 $Q$ 的三角单元。 $S_{K_Q}$ 或 $S_Q$ 和 $S_{P_0}^*$ 分别表示三角单元 $K_Q$ 和对偶单元 $K_{P_0}^*$ 的面积。设 $\mathcal{T}_h$ 和 $\mathcal{T}_h^*$ 为拟均匀的，即存在与 $h$ 无关的常数 $c_1, c_2, c_3 > 0$ ，使得

$$
c _ {1} h ^ {2} \leqslant S _ {Q} \leqslant h ^ {2}, \quad \forall Q \in \Omega_ {h} ^ {*}, \tag {7.2}
$$

![](images/13fb4eb1f46ef8d10b908dc14da5bd6c5a3ae652876ac7d721641412fb0aad07.jpg)  
图7.1 重心对偶剖分

![](images/5c2d8e40102b2f51c87529ebed2951589ff3c81fb0e100a20feebdb517949971.jpg)  
图7.2 外心对偶剖分

$$
c _ {2} h ^ {2} \leqslant S _ {P _ {0}} ^ {*} \leqslant c _ {3} h ^ {2}, \quad \forall P _ {0} \in \bar {\Omega} _ {h}. \tag {7.3}
$$

注意到条件 (7.2) 蕴含条件 (7.3). 以下总设剖分是拟均匀的.

取试探函数空间 $U_{h}$ 为相应于 $\mathcal{T}_h$ 的一次有限元空间，即

$$
U _ {h} = \left\{u _ {h}: u _ {h} \in C (\Omega), u _ {h} | _ {K} \in \mathcal {P} _ {1}, \forall K \in \mathcal {T} _ {h}, u _ {h} | _ {\partial \Omega} = 0 \right\},
$$

其中 $\mathcal{P}_1$ 为次数不超过1的多项式

对 $u \in H_0^1(\Omega)$ , 设 $\Pi_h u$ 是 $u$ 往试探函数空间 $U_h$ 的插值投影. 由 Sobolev 空间插值定理, 若 $u \in H^2(\Omega)$ , 则

$$
| u - \Pi_ {h} u | _ {m} \leqslant C h ^ {2 - m} | u | _ {2}, \quad m = 0, 1. \tag {7.4}
$$

检验函数空间 $V_{h}$ 取为相应于 $\mathcal{T}_h^*$ 的分片常数函数空间，它相应于点 $P_0\in \hat{\Omega}_h$ 的基函数为

$$
\psi_ {P _ {0}} (x, y) = \left\{ \begin{array}{l l} 1, & (x, y) \in K _ {P _ {0}} ^ {*}, \\ 0, & \text {其 他}. \end{array} \right. \tag {7.5}
$$

任一 $v_{h}\in V_{h}$ ，可唯一表示为

$$
v _ {h} = \sum_ {P _ {0} \in \hat {\Omega} _ {h}} v _ {h} \left(P _ {0}\right) \psi_ {P _ {0}}. \tag {7.6}
$$

对 $w \in U$ , 设 $\Pi_h^* w$ 是 $w$ 往检验函数空间 $V_h$ 的插值投影:

$$
\varPi_ {h} ^ {*} w = \sum_ {P _ {0} \in \hat {\Omega} _ {h}} w (P _ {0}) \psi_ {P _ {0}}. \tag {7.7}
$$

由插值理论有

$$
\left| w - \Pi_ {h} ^ {*} w \right| _ {0} \leqslant C h \left| w \right| _ {1}. \tag {7.8}
$$

# 7.1.2 线性元有限体积法

线性元有限体积格式为[19]：求 $u_{h}\in U_{h}$ ，使得

$$
a \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in V _ {h}, \tag {7.9}
$$

或等价地有

$$
a \left(u _ {h}, \psi_ {P _ {0}}\right) = (f, \psi_ {P _ {0}}), \quad \forall P _ {0} \in \mathring {\Omega} _ {h}, \tag {7.10}
$$

其中

$$
a \left(u _ {h}, v _ {h}\right) = \sum_ {P _ {0} \in \hat {\Omega} _ {h}} v _ {h} \left(P _ {0}\right) a \left(u _ {h}, \psi_ {P _ {0}}\right), \tag {7.11}
$$

$$
\begin{array}{l} a \left(u _ {h}, \psi_ {P _ {0}}\right) = - \int_ {\partial K _ {P _ {0} ^ {*}}} \left(w _ {h} ^ {(1)} n _ {x} + w _ {h} ^ {(2)} n _ {y}\right) d s \\ = - \int_ {\partial K _ {P _ {0} ^ {*}}} w _ {h} ^ {(1)} \mathrm {d} y + \int_ {\partial K _ {P _ {0} ^ {*}}} w _ {h} ^ {(2)} \mathrm {d} x. \tag {7.12} \\ \end{array}
$$

这里

$$
w _ {h} ^ {(1)} = \kappa_ {1 1} \frac {\partial u _ {h}}{\partial x} + \kappa_ {1 2} \frac {\partial u _ {h}}{\partial y}, \quad w _ {h} ^ {(2)} = \kappa_ {2 1} \frac {\partial u _ {h}}{\partial x} + \kappa_ {2 2} \frac {\partial u _ {h}}{\partial y}.
$$

上式中 $\pmb {n} = (n_x,n_y)^{\mathrm{T}}$ 为 $\partial K_{P_0^*}$ 的单位外法向量，并且

$$
n _ {x} \mathrm {d} s = \mathrm {d} y, \quad n _ {y} \mathrm {d} s = - \mathrm {d} x,
$$

另外，右端积分为

$$
(f, \psi_ {P _ {0}}) = \int_ {K _ {P _ {0}} ^ {*}} f d x d y.
$$

# 7.1.3 稳定性分析

在下面的分析中，我们假设对偶剖分始终为重心对偶剖分.首先在 $U_{h}$ 中引进离散的零模、半模和全模：

$$
\left\| u _ {h} \right\| _ {0, h} = \left(\sum_ {K \in \mathcal {T} _ {h}} \left| u _ {h} \right| _ {0, h, K} ^ {2}\right) ^ {\frac {1}{2}}, \tag {7.13}
$$

$$
\left| u _ {h} \right| _ {1, h} = \left(\sum_ {K \in \mathcal {T} _ {h}} \left| u _ {h} \right| _ {1, h, K} ^ {2}\right) ^ {\frac {1}{2}}, \tag {7.14}
$$

$$
\left\| u _ {h} \right\| _ {1, h} = \left(\left\| u _ {h} \right\| _ {0, h} ^ {2} + \left| u _ {h} \right| _ {1, h} ^ {2}\right) ^ {\frac {1}{2}}. \tag {7.15}
$$

其中 $K = K_{Q} = \triangle P_{i}P_{j}P_{k}$

$$
\begin{array}{l} \left| u _ {h} \right| _ {0, h, K} = \left[ \frac {1}{3} \left(u _ {i} ^ {2} + u _ {j} ^ {2} + u _ {k} ^ {2}\right) S _ {Q} \right] ^ {\frac {1}{2}}, \\ \left| u _ {h} \right| _ {1, h, K} = \left\{\left[ \left(\frac {\partial u _ {h} (Q)}{\partial x}\right) ^ {2} + \left(\frac {\partial u _ {h} (Q)}{\partial y}\right) ^ {2} \right] S _ {Q} \right\} ^ {\frac {1}{2}}. \\ \end{array}
$$

命题7.1 $|\cdot |_{1,h}$ 与 $|\cdot |_1$ 一致； $\| \cdot \|_{0,h}$ 和 $\| \cdot \|_{1,h}$ 分别与 $\| \cdot \| _0,\| \cdot \| _1$ 等价，即存在与 $U_{h}$ 无关的正常数 $c_{1},c_{2},c_{3},c_{4}$ ，使得

$$
c _ {1} \| u _ {h} \| _ {0, h} \leqslant \| u _ {h} \| _ {0} \leqslant c _ {2} \| u _ {h} \| _ {0, h}, \quad \forall u _ {h} \in U _ {h}, \tag {7.16a}
$$

$$
c _ {3} \| u _ {h} \| _ {1, h} \leqslant \| u _ {h} \| _ {1} \leqslant c _ {4} \| u _ {h} \| _ {1, h}, \quad \forall u _ {h} \in U _ {h}. \tag {7.16b}
$$

证明 由于 $\frac{\partial u_h}{\partial x}$ 和 $\frac{\partial u_h}{\partial y}$ 在每个单元 $K$ 内为常数，因此有

$$
\begin{array}{l} | u _ {h} | _ {1, K} ^ {2} = \iint_ {K} \nabla u _ {h} \cdot \nabla u _ {h} \mathrm {d} x \mathrm {d} y = \iint_ {K} \left[ \left(\frac {\partial u _ {h}}{\partial x}\right) ^ {2} + \left(\frac {\partial u _ {h}}{\partial y}\right) ^ {2} \right] \mathrm {d} x \mathrm {d} y \\ = \left[ \left(\frac {\partial u _ {h} (Q)}{\partial x}\right) ^ {2} + \left(\frac {\partial u _ {h} (Q)}{\partial y}\right) ^ {2} \right] S _ {Q} = | u _ {h} | _ {1, h, K} ^ {2}. \\ \end{array}
$$

即 $|\cdot |_{1,h}$ 与 $|\cdot |_1$ 一致.由于 $u_{h}$ 在 $K$ 上为一次函数，因此 $u_h^2$ 为二次函数，故利用二次精度的数值积分公式即得

$$
\begin{array}{l} \| u _ {h} \| _ {0, K} ^ {2} = \int_ {K} u _ {h} ^ {2} \mathrm {d} x \mathrm {d} y = \frac {1}{3} \left[ u _ {h} ^ {2} \left(M _ {i}\right) + u _ {h} ^ {2} \left(M _ {j}\right) + u _ {h} ^ {2} \left(M _ {k}\right) \right] S _ {Q} \\ = \frac {1}{3} \left[ \left(\frac {u _ {k} + u _ {j}}{2}\right) ^ {2} + \left(\frac {u _ {k} + u _ {i}}{2}\right) ^ {2} + \left(\frac {u _ {i} + u _ {j}}{2}\right) ^ {2} \right] S _ {Q} \\ = \frac {1}{1 2} \left[ \left(u _ {i} ^ {2} + u _ {j} ^ {2} + u _ {k} ^ {2}\right) + \left(u _ {i} + u _ {j} + u _ {k}\right) ^ {2} \right] S _ {Q}, \tag {7.17} \\ \end{array}
$$

其中 $M_{i}, M_{j}, M_{k}$ 分别为 $\overline{P_{j}P_{k}}, \overline{P_{k}P_{i}}, \overline{P_{i}P_{j}}$ 的中点，如图7.3所示

![](images/d6dec33cad3f420707d310146915fd9a8a05d6efd817db320119a4fff021a25e.jpg)  
图7.3 三角单元 $K$

由 (7.17) 有

$$
\frac {1}{1 2} \left(u _ {i} ^ {2} + u _ {j} ^ {2} + u _ {k} ^ {2}\right) S _ {Q} \leqslant \| u _ {h} \| _ {0, K} ^ {2} \leqslant \frac {1}{3} \left(u _ {i} ^ {2} + u _ {j} ^ {2} + u _ {k} ^ {2}\right) S _ {Q}.
$$

因此

$$
\frac {1}{4} | u _ {h} | _ {0, h, K} ^ {2} \leqslant \| u _ {h} \| _ {0, K} ^ {2} \leqslant | u _ {h} | _ {0, h, K} ^ {2}.
$$

对所有单元求和，即得

$$
\frac {1}{4} \| u _ {h} \| _ {0, h} ^ {2} \leqslant \| u _ {h} \| _ {0} ^ {2} \leqslant \| u _ {h} \| _ {0, h} ^ {2}.
$$

从而 $\| \cdot \|_{0,h}$ 和 $\| \cdot \|_0$ 等价, 进一步可得 $\| \cdot \|_{1,h}$ 与 $\| \cdot \|_1$ 等价. 即 (7.16) 成立.

定理7.1 当 $h$ 充分小时， $a(u_h, \Pi_h^* u_h)$ 正定，即有 $h_0 > 0, \alpha > 0$ ，使当 $0 < h \leqslant h_0$ 时，

$$
a \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) \geqslant \alpha \| u _ {h} \| _ {1} ^ {2}, \quad \forall u _ {h} \in U _ {h}. \tag {7.18}
$$

证明 由(7.7)和(7.12)有

$$
\begin{array}{l} a \left(u _ {h}, \Pi_ {h} ^ {*} \bar {u} _ {h}\right) = \sum_ {P _ {0} \in \hat {\Omega} _ {h}} \bar {u} _ {h} \left(P _ {0}\right) a \left(u _ {h}, \psi_ {P _ {0}}\right) \\ = \sum_ {P _ {0} \in \hat {\Omega} _ {h}} \bar {u} _ {h} (P _ {0}) \int_ {\partial K _ {P _ {0}} ^ {*}} \left(- w _ {h} ^ {(1)} \mathrm {d} y + w _ {h} ^ {(2)} \mathrm {d} x\right) \\ = \sum_ {K \in \mathcal {T} _ {h}} I _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} \bar {u} _ {h}\right), \tag {7.19} \\ \end{array}
$$

其中

$$
I _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} \bar {u} _ {h}\right) = \sum_ {P \in \mathring {K}} \bar {u} _ {h} (P) \int_ {\partial K _ {P} ^ {*} \cap K} \left(- w _ {h} ^ {(1)} \mathrm {d} y + w _ {h} ^ {(2)} \mathrm {d} x\right), \tag {7.20}
$$

记号 $\hat{K}$ 表示 $K = \triangle P_{i}P_{j}P_{k}$ 的三个顶点的集合

先证近似式 $a_{h}(u_{h},\Pi_{h}^{*}u_{h}) = \sum_{K\in \mathcal{T}_{h}}\tilde{I}_{K}(u_{h},\Pi_{h}^{*}u_{h})$ 的正定性，其中

$$
\begin{array}{l} \tilde {I} _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} \bar {u} _ {h}\right) = \sum_ {P \in \dot {K}} \bar {u} _ {h} (P) \int_ {\partial K _ {P} ^ {*} \cap K} (- w _ {h} ^ {(1)} (Q) d y + w _ {h} ^ {(2)} (Q) d x) \\ = \left[ w _ {h} ^ {(1)} (Q) \left(y _ {M _ {k}} - y _ {M _ {j}}\right) + w _ {h} ^ {(2)} (Q) \left(x _ {M _ {j}} - x _ {M _ {k}}\right) \right] \bar {u} _ {h} \left(P _ {i}\right) + \\ \left[ w _ {h} ^ {(1)} (Q) \left(y _ {M _ {i}} - y _ {M _ {k}}\right) + w _ {h} ^ {(2)} (Q) \left(x _ {M _ {k}} - x _ {M _ {i}}\right) \right] \bar {u} _ {h} \left(P _ {j}\right) + \\ \left[ w _ {h} ^ {(1)} (Q) \left(y _ {M _ {j}} - y _ {M _ {i}}\right) + w _ {h} ^ {(2)} (Q) \left(x _ {M _ {i}} - x _ {M _ {j}}\right) \right] \bar {u} _ {h} \left(P _ {k}\right). \\ \end{array}
$$

由(6.39)知

$$
\tilde {I} _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) = \left(w _ {h} ^ {(1)} (Q) \frac {\partial u _ {h}}{\partial x} + w _ {h} ^ {(2)} (Q) \frac {\partial u _ {h}}{\partial y}\right) S _ {Q}
$$

$$
= \left[ \kappa_ {1 1} (Q) \left(\frac {\partial u _ {h}}{\partial x}\right) ^ {2} + \left(\kappa_ {1 2} (Q) + \kappa_ {2 1} (Q)\right) \frac {\partial u _ {h}}{\partial x} \frac {\partial u _ {h}}{\partial y} + \kappa_ {2 2} (Q) \left(\frac {\partial u _ {h}}{\partial y}\right) ^ {2} \right] S _ {Q}, \tag {7.21}
$$

其中

$$
w _ {h} ^ {(1)} (Q) = \kappa_ {1 1} (Q) \frac {\partial u _ {h}}{\partial x} + \kappa_ {1 2} (Q) \frac {\partial u _ {h}}{\partial y},
$$

$$
w _ {h} ^ {(2)} (Q) = \kappa_ {2 1} (Q) \frac {\partial u _ {h}}{\partial x} + \kappa_ {2 2} (Q) \frac {\partial u _ {h}}{\partial y}.
$$

由椭圆性条件知

$$
\tilde {I} _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) \geqslant \gamma \left[ \left(\frac {\partial u _ {h} (Q)}{\partial x}\right) ^ {2} + \left(\frac {\partial u _ {h} (Q)}{\partial y}\right) ^ {2} \right] S _ {Q}.
$$

再由命题7.1并注意在 $H_0^1$ 中半模与全模的等价性知，存在 $\gamma' > 0$ 使得

$$
a _ {h} \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) \geqslant \sum_ {K \in \mathcal {T} _ {h}} \left| u _ {h} \right| _ {1, K} ^ {2} = \gamma \left| u _ {h} \right| _ {1, h} ^ {2} \geqslant \gamma^ {\prime} \| u _ {h} \| _ {1} ^ {2}, \quad \forall u _ {h} \in U _ {h}. \tag {7.22}
$$

再证 $a(u_h,\Pi_h^* u_h)$ 的正定性. 易见

$$
\begin{array}{l} I _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) - \tilde {I} _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) \\ = \sum_ {P \in \tilde {K}} \left\{- \int_ {\partial K _ {P} ^ {*} \cap K} \left[ \left(w _ {h} ^ {(1)} - w _ {h} ^ {(1)} (Q)\right) d y - \left(w _ {h} ^ {(2)} - w _ {h} ^ {(2)} (Q)\right) d x \right] \right\} u _ {h} (P) \\ = \sum_ {l = i, j, k} \int_ {\overline {{M _ {l} Q}}} \left[ \left(w _ {h} ^ {(1)} - w _ {h} ^ {(1)} (Q)\right) d y - \left(w _ {h} ^ {(2)} - w _ {h} ^ {(2)} (Q)\right) d x \right] \left(u _ {l + 2} - u _ {l + 1}\right), \tag {7.23} \\ \end{array}
$$

其中 $u_{i + 1} = u_j,u_{j + 1} = u_k,u_{k + 1} = u_i,u_l = u_h(P_l)$ .由于在 $K$ 内 $\frac{\partial u_h}{\partial x}$ 和 $\frac{\partial u_h}{\partial y}$ 为常数，故有

$$
w _ {h} ^ {(i)} - w _ {h} ^ {(i)} (Q) = \left(\kappa_ {i 1} - \kappa_ {i 1} (Q)\right) \frac {\partial u _ {h}}{\partial x} + \left(\kappa_ {i 2} - \kappa_ {i 2} (Q)\right) \frac {\partial u _ {h}}{\partial y},
$$

$$
\left| w _ {h} ^ {(i)} - w _ {h} ^ {(i)} (Q) \right| \leqslant C h \left(\left| \frac {\partial u _ {h}}{\partial x} \right| + \left| \frac {\partial u _ {h}}{\partial y} \right|\right), \quad i = 1, 2. \tag {7.24}
$$

由 Taylor 展开式及 $u_h$ 在 $K$ 内为线性函数, 从而有

$$
u _ {l + 2} - u _ {l + 1} = \frac {\partial u _ {h}}{\partial x} \left(x _ {P _ {l + 2}} - x _ {P _ {l + 1}}\right) + \frac {\partial u _ {h}}{\partial y} \left(y _ {P _ {l + 2}} - y _ {P _ {l + 1}}\right),
$$

$$
\left| u _ {l + 2} - u _ {l + 1} \right| \leqslant h \left(\left| \frac {\partial u _ {h}}{\partial x} \right| + \left| \frac {\partial u _ {h}}{\partial y} \right|\right), \quad l = i, j, k. \tag {7.25}
$$

由（7.24）（7.25）及剖分的正则性知

$$
\left| \int_ {\overline {{M _ {l} Q}}} \left[ \left(w _ {h} ^ {(1)} - w _ {h} ^ {(1)} (Q)\right) \mathrm {d} y - \left(w _ {h} ^ {(2)} - w _ {h} ^ {(2)} (Q)\right) \mathrm {d} x \right] \left(u _ {l + 2} - u _ {l + 1}\right) \right|
$$

$$
\begin{array}{l} \leqslant C h ^ {3} \left(\left| \frac {\partial u _ {h}}{\partial x} \right| + \left| \frac {\partial u _ {h}}{\partial y} \right|\right) ^ {2} \\ \leqslant \tilde {C} h \left[ \left(\frac {\partial u _ {h}}{\partial x}\right) ^ {2} + \left(\frac {\partial u _ {h}}{\partial y}\right) ^ {2} \right] S _ {Q} \\ = \tilde {C} h \left| u _ {h} \right| _ {1, K} ^ {2}. \tag {7.26} \\ \end{array}
$$

由（7.23）（7.26）和命题7.1知

$$
\begin{array}{l} \left| a \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) - a _ {h} \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) \right| \\ = \left| \sum_ {K \in \mathcal {T} _ {h}} \left[ I _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) - \tilde {I} _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) \right] \right| \leqslant \hat {C h} \| u _ {h} \| _ {1} ^ {2}. \tag {7.27} \\ \end{array}
$$

联立 (7.22) 和 (7.27) 即得 (7.18).

# 7.1.4 误差估计

定义在 $\mathcal{T}_h$ 上的分片 $H^2$ 空间

$$
H _ {h} ^ {2} (\Omega) = \{u \in C (\Omega): u | _ {K} \in H ^ {2} (K), \forall K \in \mathcal {T} _ {h} \}
$$

及其上的范数

$$
| u | _ {2, h} := \left(\sum_ {K \in \mathcal {T} _ {h}} | u | _ {2, K} ^ {2}\right) ^ {\frac {1}{2}}, \quad \forall u \in H _ {h} ^ {2} (\Omega).
$$

引理7.1 若 $u \in H_h^2(\Omega)$ 且 $u_h \in U_h$ ，则存在一个正常数 $C$ 使得

$$
\left| a _ {h} \left(u, \Pi_ {h} ^ {*} u _ {h}\right) \right| \leqslant C \left(\left| u \right| _ {1} + h \left| u \right| _ {2, h}\right) \left| u _ {h} \right| _ {1}. \tag {7.28}
$$

证明 对于每个单元 $K \in \mathcal{T}_h$ ，记在 $K$ 内所有对偶单元的边界线段集合为 $L_K^*$ ，如图7.3所示， $L_K^* = \{\overline{QM_i}, \overline{QM_j}, \overline{QM_k}\}$ 。因此，可以整理得

$$
\begin{array}{l} a _ {h} (u, \Pi_ {h} ^ {*} u _ {h}) = - \sum_ {K \in \mathcal {T} _ {h}} \sum_ {K _ {P} ^ {*} \in \mathcal {T} _ {h} ^ {*}} \int_ {\partial K _ {P} ^ {*} \cap K} (\kappa \nabla u) \cdot \boldsymbol {n} (\Pi_ {h} ^ {*} u _ {h}) (P) d s \\ = - \sum_ {K \in \mathcal {T} _ {h}} \sum_ {l ^ {*} \in L _ {K} ^ {*}} \int_ {l ^ {*}} (\kappa \nabla u) \cdot \boldsymbol {n} [ \Pi_ {h} ^ {*} u _ {h} ] _ {l ^ {*}} d s, \\ \end{array}
$$

其中 $\pmb{n}$ 是在线段 $l^{*}$ 处的单位外法向量，由一个对偶单元 $K_{1}^{*}$ 指向另一个对偶单元 $K_{2}^{*}$ 并且

$$
\left[ \Pi_ {h} ^ {*} u _ {h} \right] _ {l ^ {*}} := \Pi_ {h} ^ {*} u _ {h} | _ {K _ {1} ^ {*}} - \Pi_ {h} ^ {*} u _ {h} | _ {K _ {2} ^ {*}}.
$$

例如, 若 $l^{*} = \overline{QM_{k}}$ 处单位外法向量 $\pmb{n}$ 由 $K_{P_i}^*$ 指向 $K_{P_j}^*$ , 则 $[\varPi_h^* u_h]_{l^*} = u_h(P_i) - u_h(P_j)$ .

由Cauchy-Schwarz不等式，可得

$$
| a _ {h} (u, \Pi_ {h} ^ {*} u _ {h}) | ^ {2} \leqslant C \Bigg (\sum_ {K \in \mathcal {T} _ {h}} \sum_ {l ^ {*} \in L _ {K} ^ {*}} ([ \Pi_ {h} ^ {*} u _ {h} ] _ {l ^ {*}}) ^ {2} \Bigg) \Bigg (\sum_ {K \in \mathcal {T} _ {h}} \sum_ {l ^ {*} \in L _ {K} ^ {*}} | l ^ {*} | \int_ {l ^ {*}} (\kappa \nabla u \cdot \pmb {n}) ^ {2} \mathrm {d} s \Bigg). \tag {7.29}
$$

一方面，由（7.25）及剖分的正则性有

$$
\begin{array}{l} \sum_ {K \in \mathcal {T} _ {h}} \sum_ {l ^ {*} \in L _ {K} ^ {*}} \left(\left[ \Pi_ {h} ^ {*} u _ {h} \right] _ {l ^ {*}}\right) ^ {2} = \sum_ {K \in \mathcal {T} _ {h}} \left[ \left(u _ {h} \left(P _ {i}\right) - u _ {h} \left(P _ {j}\right)\right) ^ {2} + \left(u _ {h} \left(P _ {i}\right) - u _ {h} \left(P _ {k}\right)\right) ^ {2} + \left(u _ {h} \left(P _ {j}\right) - u _ {h} \left(P _ {k}\right)\right) ^ {2} \right] \\ \leqslant C \sum_ {K \in \mathcal {T} _ {h}} h ^ {2} \left(\left| \frac {\partial u _ {h}}{\partial x} \right| ^ {2} + \left| \frac {\partial u _ {h}}{\partial y} \right| ^ {2}\right) \\ \leqslant C \left| u _ {h} \right| _ {1} ^ {2}. \tag {7.30} \\ \end{array}
$$

另一方面，由于 $\kappa$ 是对称矩阵且每个元素充分光滑，可以得到

$$
\sum_ {K \in \mathcal {T} _ {h}} \sum_ {l ^ {*} \in L _ {K} ^ {*}} | l ^ {*} | \int_ {l ^ {*}} ((\kappa \nabla u) \cdot n) ^ {2} \mathrm {d} s \leqslant C h \sum_ {K \in \mathcal {T} _ {h}} \sum_ {l ^ {*} \in L _ {K} ^ {*}} \int_ {l ^ {*}} | \nabla u | ^ {2} \mathrm {d} s. \tag {7.31}
$$

记 $\varphi = \nabla u$ ，由单元 $K$ 与参考单元 $\hat{K}$ 之间的关系可得

$$
\sum_ {l ^ {*} \in L _ {K} ^ {*}} \int_ {l ^ {*}} | \varphi | ^ {2} \mathrm {d} s \leqslant C h \sum_ {\hat {l} ^ {*} \in L _ {\hat {K}} ^ {*}} \int_ {\hat {l} ^ {*}} | \hat {\varphi} | ^ {2} \mathrm {d} \hat {s}. \tag {7.32}
$$

根据迹定理，有

$$
\sum_ {\hat {l} ^ {*} \in L _ {\hat {K}} ^ {*}} \int_ {\hat {l} ^ {*}} | \hat {\varphi} | ^ {2} \mathrm {d} \hat {s} \leqslant C \| \hat {\varphi} \| _ {1, \hat {K}} ^ {2}, \tag {7.33}
$$

注意到

$$
\| \hat {\varphi} \| _ {0, \hat {K}} ^ {2} \leqslant C h ^ {- 2} \| \varphi \| _ {0, K} ^ {2}, \quad | \hat {\varphi} | _ {1, \hat {K}} ^ {2} \leqslant C | \varphi | _ {1, K} ^ {2}. \tag {7.34}
$$

结合 (7.32) (7.33) 和 (7.34) 可得

$$
\sum_ {l ^ {*} \in L _ {K} ^ {*}} \int_ {l ^ {*}} | \varphi | ^ {2} d s \leqslant C (h ^ {- 1} \| \varphi \| _ {0, K} ^ {2} + h | \varphi | _ {1, K} ^ {2}),
$$

将上式带入（7.31）得到估计式

$$
\sum_ {K \in \mathcal {T} _ {h}} \sum_ {l ^ {*} \in L _ {K} ^ {*}} | l ^ {*} | \int_ {l ^ {*}} (\kappa \nabla u \cdot n) ^ {2} d s \leqslant C (| u | _ {1} ^ {2} + h ^ {2} | u | _ {2, h} ^ {2}). \tag {7.35}
$$

最后，由（7.29）和（7.35）可直接推出连续性结果(7.28)，引理得证

定理7.2 设 $u$ 是问题（7.1）的广义解， $u_h$ 是有限体积元格式（7.10）的解。若 $u \in H^2(\Omega)$ ，则有误差估计：

$$
\left\| u - u _ {h} \right\| _ {1} \leqslant C h | u | _ {2},
$$

$$
\| u - u _ {h} \| _ {0} \leqslant C h ^ {2} | u | _ {3}. \tag {7.36}
$$

证明 显然有

$$
a \left(u - u _ {h}, \psi_ {P _ {0}}\right) = 0, \quad \forall P _ {0} \in \hat {\Omega} _ {h}, \tag {7.37}
$$

由定理7.1和（7.37）有

$$
\begin{array}{l} \left\| u _ {h} - \Pi_ {h} u \right\| _ {1} ^ {2} \leqslant \frac {1}{\alpha} a \left(u _ {h} - \Pi_ {h} u, \Pi_ {h} ^ {*} \left(u _ {h} - \Pi_ {h} u\right)\right) \\ = \frac {1}{\alpha} a \left(u _ {h} - u + u - \Pi_ {h} u, \Pi_ {h} ^ {*} \left(u _ {h} - \Pi_ {h} u\right)\right) \\ = \frac {1}{\alpha} a (u - \Pi_ {h} u, \Pi_ {h} ^ {*} (u _ {h} - \Pi_ {h} u)), \\ \end{array}
$$

从而由引理7.1可得

$$
\begin{array}{l} \left\| u _ {h} - \Pi_ {h} u \right\| _ {1} ^ {2} \leqslant \frac {1}{\alpha} \left| a \left(u - \Pi_ {h} u, \Pi_ {h} ^ {*} \left(u _ {h} - \Pi_ {h} u\right)\right) \right| \\ \leqslant C \left(| u - \Pi_ {h} u | _ {1} + h | u | _ {2}\right) | u _ {h} - \Pi_ {h} u | _ {1}. \tag {7.38} \\ \end{array}
$$

利用 (7.4) 得

$$
| u - \Pi_ {h} u | _ {1} \leqslant C h | u | _ {2}, \tag {7.39}
$$

联立 (7.38) 和 (7.39) 可推出

$$
\left\| u _ {h} - \Pi_ {h} u \right\| _ {1} \leqslant C h | u | _ {2}.
$$

其中 $H^{1}(I)$ 空间里半模和全模是等价的.再由三角不等式知

$$
\| u - u _ {h} \| _ {1} \leqslant \| u - \Pi_ {h} u \| _ {1} + \| u _ {h} - \Pi_ {h} u \| _ {1}.
$$

利用（7.39）知（7.2）成立，从而定理得证

注7.1 若 $\kappa$ 为常数矩阵，根据 $I_K(u_h, \Pi_h^* \bar{u}_h)$ 和 $\tilde{I}_K(u_h, \Pi_h^* \bar{u}_h)$ 的定义，我们有 $I_K = \tilde{I}_K$ ，于是

$$
\begin{array}{l} \sum_ {K \in \mathcal {T} _ {h}} I _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} \bar {u} _ {h}\right) = \sum_ {K \in \mathcal {T} _ {h}} \tilde {I} _ {K} \left(u _ {h}, \Pi_ {h} ^ {*} \bar {u} _ {h}\right) \\ = \sum_ {K \in \mathcal {T} _ {h}} \left(w _ {h} ^ {(1)} (Q) \frac {\partial \bar {u} _ {h}}{\partial x} (Q) + w _ {h} ^ {(2)} (Q) \frac {\partial \bar {u} _ {h}}{\partial y} (Q)\right) S _ {Q} \\ = \sum_ {K \in \mathcal {T} _ {h}} \iint_ {K} (\kappa (Q) \nabla u _ {h}) \cdot \nabla \bar {u} _ {h} d x d y \\ = \iint_ {\Omega} (\kappa \nabla u _ {h}) \cdot \nabla \bar {u} _ {h} d x d y \\ = a \left(u _ {h}, \bar {u} _ {h}\right). \\ \end{array}
$$

从而成立

$$
a \left(u _ {h}, \Pi_ {h} ^ {*} \bar {u} _ {h}\right) = a \left(\bar {u} _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) = a \left(u _ {h}, \bar {u} _ {h}\right),
$$

其中 $a(u_h,\bar{u}_h)$ 是有限元法的双线性形式.此时，有限体积元法双线性形式对称，且与有限元法的双线性形式相等.

注7.2 若 $\pmb{A}$ 为变系数矩阵，有限体积元法的双线性形式不对称，即 $a(u_h, \Pi_h^* \bar{u}_h) \neq a(\bar{u}_h, \Pi_h^* u_h)$ ，但我们可以用如下方法将其对称化。

(i) 在重心对偶剖分情形下, 用变系数 $a_{ij}$ 在重心处的值 $a_{ij}(Q)$ 作数值积分处理;  
(ii) 用变系数在三角形单元上的积分平均值, 将变系数常数化, 即取

$$
\bar {a} _ {i j} := \frac {1}{| K |} \iint_ {K} a _ {i j} \mathrm {d} x \mathrm {d} y, \quad i, j = 1, 2.
$$

# 7.2 四边形网格上的有限体积元法

# 7.2.1 四边形网格剖分及对偶剖分

将 $\bar{\Omega}$ 分割成有限个严格凸的四边形之和，使不同的四边形无公共的内点，任一四边形的顶点不属于其他四边形边的内部，且边界的任一角点都是某一四边形的顶点。每个四边形称为单元，记为 $K$ 。所有单元构成 $\bar{\Omega}$ 的四边形剖分，记为 $T_{h}$ ， $h$ 表示所有四边形的最大直径。四边形的顶点称为剖分的节点。另外，当求解区域形状规则简单时，也常常对区域作矩形网格剖分，在矩形网格上构造的有限体积格式更简单。

再作和 $T_{h}$ 相应的对偶剖分.如图7.4,设 $P_0$ 是剖分 $T_{h}$ 的任一节点， $P_{i}$ $(i = 1,2,3,4)$

![](images/d1d4f3e13aa03ee1675f181f82095b87a55fa4c2ae8bb276352498920ebc286f.jpg)  
图7.4 四边形单元 $K$ 的对偶单元

是与 $P_{0}$ 相邻的节点, $M_{i}$ 是 $\overline{P_0P_i}$ 的中点, $P_{i,i + 1}(P_{45} = P_{41} = P_{14})$ 是以 $\overline{P_0P_i}$ 和 $\overline{P_0P_{i + 1}}$ 为边的四边形中与 $P_{0}$ 相对的顶点. 在四边形 $P_{0}P_{i}P_{i,i + 1}P_{i + 1}$ 内取平均中心 $Q_{i}$ ( $i = 1,2,3,4$ ), 即对边中点连线的交点. 依次连接 $M_{1}, Q_{1}, M_{2}, Q_{2}, \dots, M_{4}, Q_{4}, M_{1}$ , 得到一个围绕 $P_{0}$ 的多边形域 $K_{P_0}^*$ , 称为对偶单元. 所有对偶单元构成 $\bar{\Omega}$ 的一个新剖分 $T_h^*$ , 称为对偶剖分. $\bar{\Omega}_h$ 表示 $T_h$ 的节点集合; $\mathring{\Omega}_h = \bar{\Omega}_h \setminus \partial \Omega$ 表示内节点集合; $\Omega_h^*$ 表示对偶剖分 $T_h^*$ 的节点集合. 对于 $Q \in \Omega_h^*$ , 以 $K_Q$ 表示以 $Q$ 为平均中心的四边形单元, $S_Q$ (或 $S_{K_Q}$ ) 和 $S_P^*$ 分别表示 $K_Q$ 和 $K_P^*$ 的面积. 总设剖分 $T_h$ 和 $T_h^*$ 是拟均匀剖分, 即存在与 $h$ 无关的正常数 $c_1, c_2, c_3$ , 使得

$$
c _ {1} h ^ {2} \leqslant S _ {Q} \leqslant h ^ {2}, \quad Q \in \Omega_ {h} ^ {*}, \tag {7.40a}
$$

$$
c _ {2} h ^ {2} \leqslant S _ {P _ {0}} ^ {*} \leqslant c _ {3} h ^ {2}, \quad P _ {0} \in \bar {\Omega} _ {h}. \tag {7.40b}
$$

取 $(\xi, \eta)$ 平面上的单位正方形 $\hat{K} = \hat{E} = [0,1] \times [0,1]$ 作为参考单元. 对于任意凸四边形单元 $K_{Q} = \square P_{1}P_{2}P_{3}P_{4}, P_{i} = (x_{i},y_{i}) (i = 1,2,3,4)$ . 存在唯一的可逆双线性变换:

$$
F _ {K _ {Q}}: \left\{ \begin{array}{l} x = x _ {1} + a _ {1} \xi + a _ {2} \eta + a _ {3} \xi \eta , \\ y = y _ {1} + b _ {1} \xi + b _ {2} \eta + b _ {3} \xi \eta , \end{array} \right. \tag {7.41}
$$

其中

$$
a _ {1} = x _ {2} - x _ {1}, \quad a _ {2} = x _ {3} - x _ {1}, \quad a _ {3} = x _ {4} - x _ {3} - x _ {2} + x _ {1};
$$

$$
b _ {1} = y _ {2} - y _ {1}, \quad b _ {2} = y _ {3} - y _ {1}, \quad b _ {3} = y _ {4} - y _ {3} - y _ {2} + y _ {1}.
$$

将 $\hat{K}$ 变成 $K_{Q}$ (参见图7.5). 用 $\mathcal{J}_K$ 表示等参双线性变换 $F_{K_Q}$ 的Jacobi矩阵, 其行列式用 $J_{K}$ 表示, 则

$$
\mathcal {J} _ {K} = \left[ \begin{array}{l l} \frac {\partial x}{\partial \xi} & \frac {\partial x}{\partial \eta} \\ \frac {\partial y}{\partial \xi} & \frac {\partial y}{\partial \eta} \end{array} \right] = \left[ \begin{array}{l l} a _ {1} + a _ {3} \eta & a _ {2} + a _ {3} \xi \\ b _ {1} + b _ {3} \eta & b _ {2} + b _ {3} \xi \end{array} \right].
$$

由反函数的微分法则，得

![](images/0d58b6ac5b9175ce11b620ebf991ba980f5e257186c73ab2177520a4bc383688.jpg)  
图7.5 参考单元 $\hat{K}$ 与四边形单元 $K_{Q}$

$$
\left\{ \begin{array}{l l} \frac {\partial \xi}{\partial x} = \frac {b _ {2} + b _ {3} \xi}{J _ {K}}, & \frac {\partial \xi}{\partial y} = - \frac {a _ {2} + a _ {3} \xi}{J _ {K}}, \\ \frac {\partial \eta}{\partial x} = - \frac {b _ {1} + b _ {3} \eta}{J _ {K}}, & \frac {\partial \eta}{\partial y} = \frac {a _ {1} + a _ {3} \eta}{J _ {K}}. \end{array} \right.
$$

注7.3 变换 $F_{K_Q}$ 是一一对应的，亦即它的Jacobi行列式处处不为零：

$$
\begin{array}{l} J _ {K} (\xi , \eta) = \left| \begin{array}{l l} \frac {\partial x}{\partial \xi} & \frac {\partial x}{\partial \eta} \\ \frac {\partial y}{\partial \xi} & \frac {\partial y}{\partial \eta} \end{array} \right| = \left| \begin{array}{l l} a _ {1} + a _ {3} \eta & a _ {2} + a _ {3} \xi \\ b _ {1} + b _ {3} \eta & b _ {2} + b _ {3} \xi \end{array} \right| \\ = \left(a _ {1} b _ {2} - a _ {2} b _ {1}\right) + \left(a _ {1} b _ {3} - a _ {3} b _ {1}\right) \xi + \left(a _ {3} b _ {2} - a _ {2} b _ {3}\right) \eta , \\ \end{array}
$$

它是 $\xi, \eta$ 的一次函数，要它不为零，即在四个顶点处 $(0,0),(1,0),(0,1)$ ， $(1,1)$ 的值必须同号：

$$
J (0, 0) = a _ {1} b _ {2} - a _ {2} b _ {1} = \left| \overline {{P _ {1} P _ {2}}} \right| \left| \overline {{P _ {1} P _ {3}}} \right| \sin \angle P _ {3} P _ {1} P _ {2},
$$

即当 $0 < \angle P_{3}P_{1}P_{2} < \pi$ 时, $J(0,0) > 0$ . 在其他三点处可得到类似的结果, 亦即当 $K$ 为凸四边形时 $\Leftrightarrow J(\xi, \eta) > 0$ . 这也是为什么要求剖分 $T_{h}$ 的四边形为凸四边形.

注7.4 当 $K_{Q}$ 是一个平行四边形（包括矩形）时，我们有 $a_3 = b_3 = 0$ ，变换 $F_{K_Q}$ 为线性变换.

# 7.2.2 试探函数空间和检验函数空间

定义参考单元 $\hat{K}$ 上的双线性函数

$$
P _ {\hat {K}} (\xi , \eta) = c _ {0} + c _ {1} \xi + c _ {2} \eta + c _ {3} \xi \eta ,
$$

其中

$$
c _ {0} = u _ {P _ {1}}, \quad c _ {1} = u _ {P _ {2}} - u _ {P _ {1}}, \quad c _ {2} = u _ {P _ {3}} - u _ {P _ {1}},
$$

$$
c _ {3} = u _ {P _ {4}} - u _ {P _ {3}} - u _ {P _ {2}} + u _ {P _ {1}}.
$$

定义试探函数空间

$$
U _ {h} = \left\{u _ {h} \in C ^ {0} (\bar {\Omega}): u _ {h} | _ {K} = P _ {\hat {K}} \circ F _ {K} ^ {- 1}, \quad K \in T _ {h}, \quad P _ {\hat {K}} \in P _ {1, 1}, \quad u _ {h} | _ {\partial \Omega} = 0 \right\},
$$

其中 $P_{1,1}$ 是 $\hat{K}$ 上双线性函数构成的集合. 对任一 $u_h \in U_h$ , 在 $K_Q$ 上, 有

$$
\begin{array}{l} u _ {h} = u _ {1} (1 - \xi) (1 - \eta) + u _ {2} \xi (1 - \eta) + u _ {3} (1 - \xi) \eta + u _ {4} \xi \eta \\ = u _ {1} + \left(u _ {2} - u _ {1}\right) \xi + \left(u _ {3} - u _ {1}\right) \eta + \left(u _ {4} - u _ {3} - u _ {2} + u _ {1}\right) \xi \eta . \tag {7.42} \\ \end{array}
$$

检验函数空间 $V_{h}$ 取为相应于对偶剖分 $T_{h}^{*}$ 的分片常数空间, 其基函数如下: 对 $P_0 \in \mathring{\partial}_h$

$$
\psi_ {P _ {0}} (P) = \left\{ \begin{array}{l l} 1, & P \in K _ {P _ {0}} ^ {*}, \\ 0, & P \notin K _ {P _ {0}} ^ {*}. \end{array} \right.
$$

对任一 $v_{h}\in V_{h}$ ，有

$$
v _ {h} = \sum_ {P _ {0} \in \dot {\Omega} _ {h}} v _ {h} (P _ {0}) \psi_ {P _ {0}}.
$$

# 7.2.3 等参双线性有限体积元法

取前述的试探函数空间和检验函数空间，相应于问题（7.1）的有限体积元法为：求 $u_{h}\in U_{h}$ ，使得

$$
a \left(u _ {h}, \psi_ {P _ {0}}\right) = \left(f, \psi_ {P _ {0}}\right), \quad \forall P _ {0} \in \mathring {\Omega} _ {h}, \tag {7.43}
$$

其中

$$
\begin{array}{l} a \left(u _ {h}, v _ {h}\right) = \sum_ {P _ {0} \in \mathring {\Omega} _ {h}} v _ {h} \left(P _ {0}\right) a \left(u _ {h}, \psi_ {P _ {0}}\right), (7.44) \\ a \left(u _ {h}, \psi_ {P _ {0}}\right) = - \int_ {\partial K _ {P _ {0} ^ {*}}} \left(w _ {h} ^ {(1)} n _ {x} + w _ {h} ^ {(2)} n _ {y}\right) d s \\ = - \int_ {\partial K _ {P _ {0} ^ {*}}} w _ {h} ^ {(1)} \mathrm {d} y + \int_ {\partial K _ {P _ {0} ^ {*}}} w _ {h} ^ {(2)} \mathrm {d} x. (7.45) \\ \end{array}
$$

注7.5 在推导格式时，单元 $K$ （见图7.6）上的线积分计算需要用到如下公式：

$$
\frac {\partial u _ {h}}{\partial x} = \frac {\partial u _ {h}}{\partial \xi} \frac {\partial \xi}{\partial x} + \frac {\partial u _ {h}}{\partial \eta} \frac {\partial \eta}{\partial x},
$$

$$
\frac {\partial u _ {h}}{\partial y} = \frac {\partial u _ {h}}{\partial \xi} \frac {\partial \xi}{\partial y} + \frac {\partial u _ {h}}{\partial \eta} \frac {\partial \eta}{\partial y},
$$

在 $\overline{M_1M_3}$ 上，

$$
\mathrm {d} x = \mathrm {d} \left(x _ {1} + \frac {a _ {1}}{2} + a _ {2} \eta + \frac {a _ {3}}{2} \eta\right) = \left(a _ {2} + \frac {a _ {3}}{2}\right) \mathrm {d} \eta ,
$$

$$
\mathrm {d} y = \mathrm {d} \left(y _ {1} + \frac {b _ {1}}{2} + b _ {2} \eta + \frac {b _ {3}}{2} \eta\right) = \left(b _ {2} + \frac {b _ {3}}{2}\right) \mathrm {d} \eta ,
$$

在 $\overline{M_2M_4}$ 上，

$$
\mathrm {d} x = \mathrm {d} \left(x _ {1} + a _ {1} \xi + \frac {a _ {2}}{2} + \frac {a _ {3}}{2} \xi\right) = \left(a _ {1} + \frac {a _ {3}}{2}\right) \mathrm {d} \xi ,
$$

$$
\mathrm {d} y = \mathrm {d} \left(y _ {1} + b _ {1} \xi + \frac {b _ {2}}{2} + \frac {b _ {3}}{2} \xi\right) = \left(b _ {1} + \frac {b _ {3}}{2}\right) \mathrm {d} \xi .
$$

![](images/7301b9b171d598dc1c212df10c8436e651d34396b4d53a1695a2b0e0bb4b320b.jpg)  
图7.6 四边形单元 $K_{Q}$ 上的控制体积

# 7.2.4 收敛性

在 $U_{h}$ 上定义离散半模：

$$
| u _ {h} | _ {1, h} = \left(\sum_ {Q \in \Omega_ {h} ^ {*}} | u _ {h} | _ {1, K _ {Q}, h} ^ {2}\right) ^ {\frac {1}{2}}, \quad \forall u _ {h} \in U _ {h},
$$

其中

$$
\left| u _ {h} \right| _ {1, K _ {Q}, h} ^ {2} = \left(u _ {P _ {2}} - u _ {P _ {1}}\right) ^ {2} + \left(u _ {P _ {4}} - u _ {P _ {2}}\right) ^ {2} + \left(u _ {P _ {4}} - u _ {P _ {3}}\right) ^ {2} + \left(u _ {P _ {3}} - u _ {P _ {1}}\right) ^ {2},
$$

$$
u _ {P _ {i}} = u _ {h} (P _ {i}), i = 1, 2, 3, 4.
$$

命题7.2 半模 $|u_h|_{1,h}$ 与 $|u_h|_1$ 等价，即有常数 $\beta_{1},\beta_{2} > 0$ ，使得

$$
\beta_ {1} \left| u _ {h} \right| _ {1, h} \leqslant \left| u _ {h} \right| _ {1} \leqslant \beta_ {2} \left| u _ {h} \right| _ {1, h}.
$$

定理7.3 当 $h$ 充分小时, $a(u_h,\Pi_h^* u_h)$ 正定, 即有 $h_0 > 0, \alpha > 0$ , 使当 $0 < h \leqslant h_0$ 时,

$$
a \left(u _ {h}, \Pi_ {h} ^ {*} u _ {h}\right) \geqslant \alpha \| u _ {h} \| _ {1} ^ {2}, \quad \forall u _ {h} \in U _ {h}. \tag {7.46}
$$

由稳定性(7.46)立即可得有限体积法(7.43)有唯一解

定理7.4 设 $u \in H_0^1(\Omega) \cap H^2(\Omega)$ 为问题(7.1)的广义解, $u_h \in U_h$ 为有限体积格式（7.43）的解，且剖分条件成立，则按 $H^1$ 模有误差估计

$$
\| u - u _ {h} \| _ {1} \leqslant C h \| u \| _ {2}. \tag {7.47}
$$

若还有 $u \in H^3(\Omega)$ , 则按 $L^2$ 模有误差估计

$$
\| u - u _ {h} \| _ {0} \leqslant C h ^ {2} (\| u \| _ {2} + \| f \| _ {1}). \tag {7.48}
$$

这些理论结果的证明可以参看[20, 21], 此处略去.

# 7.2.5 数值算例

例7.1 考虑椭圆型方程

$$
\left\{ \begin{array}{l l} - \nabla (\kappa (x, y) \nabla u) u = f (x, y), & (x, y) \in \Omega = (0, 1) \times (0, 1), \\ u = 0, & (x, y) \in \Gamma = \partial \Omega , \end{array} \right.
$$

其中

$$
\kappa (x, y) = \left[ \begin{array}{c c} 1 + \mathrm {e} ^ {x + y} & 0 \\ 0 & 1 + \mathrm {e} ^ {x + y} \end{array} \right],
$$

选取右端 $f$ 使得其精确解为 $u(x,y) = \sin (2\pi x)\sin (3\pi y).$

对区域 $\Omega = (0,1)\times (0,1)$ 作正方形网格剖分， $h = \frac{1}{4},\frac{1}{8},\frac{1}{16},\frac{1}{32},\frac{1}{64},\frac{1}{128}$ 表示正方形单元的边长，用双线性元有限体积法求解上述问题.数值结果如表7.1,其中 $x_{1} = \frac{1}{4}$ $x_{2} = \frac{1}{2},x_{3} = \frac{3}{4},y_{1} = \frac{1}{4},y_{1} = \frac{1}{4},y_{2} = \frac{1}{2},y_{3} = \frac{3}{4},$ 取四位小数计算

表 7.1 例 7.1 的数值结果  

<table><tr><td>h</td><td>(xj,yk)</td><td>x1</td><td>x2</td><td>x3</td></tr><tr><td>1/4</td><td>y1</td><td>1.0728</td><td>-0.0040</td><td>-1.0719</td></tr><tr><td>1/8</td><td></td><td>0.7941</td><td>0.0003</td><td>-0.7939</td></tr><tr><td>1/16</td><td></td><td>0.7283</td><td>0.0001</td><td>-0.7282</td></tr><tr><td>1/32</td><td></td><td>0.7124</td><td>0.0000</td><td>-0.7124</td></tr><tr><td>1/64</td><td></td><td>0.7084</td><td>0.0000</td><td>-0.7084</td></tr><tr><td>1/128</td><td></td><td>0.7074</td><td>0.0000</td><td>-0.7074</td></tr><tr><td>精确解</td><td></td><td>0.7071</td><td>0.0000</td><td>-0.7071</td></tr><tr><td>1/4</td><td>y2</td><td>-1.5268</td><td>0.0052</td><td>1.5273</td></tr><tr><td>1/8</td><td></td><td>-1.1238</td><td>-0.0005</td><td>1.1237</td></tr><tr><td>1/16</td><td></td><td>-1.0301</td><td>-0.0002</td><td>1.0301</td></tr><tr><td>1/32</td><td></td><td>-1.0075</td><td>-0.0000</td><td>1.0075</td></tr><tr><td>1/64</td><td></td><td>-1.0019</td><td>-0.0000</td><td>1.0019</td></tr><tr><td>1/128</td><td></td><td>-1.0005</td><td>-0.0000</td><td>1.0005</td></tr><tr><td>精确解</td><td></td><td>-1.0000</td><td>0.0000</td><td>1.0000</td></tr><tr><td>1/4</td><td>y3</td><td>1.0886</td><td>-0.0032</td><td>-1.0898</td></tr><tr><td>1/8</td><td></td><td>0.7954</td><td>0.0005</td><td>-0.7954</td></tr><tr><td>1/16</td><td></td><td>0.7285</td><td>0.0002</td><td>-0.7285</td></tr><tr><td>1/32</td><td></td><td>0.7124</td><td>0.0000</td><td>-0.7124</td></tr><tr><td>1/64</td><td></td><td>0.7084</td><td>0.0000</td><td>-0.7084</td></tr><tr><td>1/128</td><td></td><td>0.7074</td><td>0.0000</td><td>-0.7074</td></tr><tr><td>精确解</td><td></td><td>0.7071</td><td>0.0000</td><td>-0.7071</td></tr></table>

# 7.2.6 习题

考虑Possion方程边值问题

$$
\left\{ \begin{array}{l l} - \Delta u = f, & (x, y) \in \Omega , \\ u | _ {\partial \Omega} = 0, & (x, y) \in \Gamma = \partial \Omega , \end{array} \right.
$$

其中 $\Omega$ 为平面矩形区域, $\Gamma = \partial \Omega$ 为其边界, 且 $f \in L^{2}(\Omega)$ . 对求解区域 $\Omega$ 作矩形网格剖分, 通过单元分析法, 估计双线性有限元法与双线性有限体积法二者双线性形式的差.

![](images/a96c6f753bec9b702ba44324acd4b1325b60f03634c82af02ac028944b247729.jpg)

# 第8章

# 间断 Galerkin 法

同第6章介绍的经典有限元方法相比，本章介绍的间断有限元方法[36, 37, 38, 39]更加灵活，不要求逼近空间中的函数跨过相邻单元时保持连续性，可以用于更一般类型的网格剖分，允许同一网格中的不同单元采用不同种类的逼近函数空间；具有局部守恒性质；但自由度数往往相对较多。间断有限元方法同样可以应用于非常多种类的偏微分方程。本章以如下模型问题为例来介绍几种间断有限元法：

$$
\left\{ \begin{array}{l l} - \nabla \cdot (\kappa (\boldsymbol {x}) \nabla u) = f (\boldsymbol {x}), & \boldsymbol {x} \in \Omega , \\ u = 0, & \boldsymbol {x} \in \Gamma_ {D}, \\ \kappa \nabla u \cdot \boldsymbol {n} = 0, & \boldsymbol {x} \in \Gamma_ {N}, \end{array} \right. \tag {8.1}
$$

其中 $\Omega \subset \mathbb{R}^d$ 是多边形 (或多面体) 区域, $\Gamma \coloneqq \partial \Omega = \Gamma_D \cup \Gamma_N$ . 假设 $f \in L^{2}(\Omega)$ , 则存在常数 $a_1 \geqslant a_0 > 0$ 使得矩阵 $\kappa$ 满足

$$
a _ {0} (\boldsymbol {\xi}, \boldsymbol {\xi}) \leqslant (\kappa (\boldsymbol {x}) \boldsymbol {\xi}, \boldsymbol {\xi}) \leqslant a _ {1} (\boldsymbol {\xi}, \boldsymbol {\xi}), \quad \forall \boldsymbol {\xi} \in \mathbb {R} ^ {d}, \boldsymbol {x} \in \Omega .
$$

# 8.1 内罚间断 Galerkin 法

内罚间断 Galerkin (IPDG) 法 [36, 39] 是对数值解的梯度在单元公共边 (面) 上的跳量加罚, 在对流扩散问题及声波散射问题等领域中有应用.

# 8.1.1 离散格式

设 $\mathcal{M}_h$ 是 $\Omega$ 的一个三角剖分. 在本章中, 为了叙述简单, “单元的边”在二维情形就是本来的含义, 在三维情形指四面体单元的面. $\mathcal{E}_h$ 为所有单元边的集合. $\mathcal{E}_h^I$ 是 $\Omega$ 内部的单元边的集合, $\mathcal{E}_h^D$ 是 $\Gamma_D$ 上的单元边的集合, $\mathcal{E}_h^N$ 是 $\Gamma_N$ 上的单元边的集合. 记 $\mathcal{E}_h^{ID} = \mathcal{E}_h^I \cup \mathcal{E}_h^D$ , $\mathcal{E}_h^{IN} = \mathcal{E}_h^I \cup \mathcal{E}_h^N$ . 对任意 $K \in \mathcal{M}_h$ , 记 $h_K$ 为 $K$ 的直径. $\forall e \in \mathcal{E}_h$ , 记 $h_e$ 为边 $e$ 的直径. 记 $(\cdot, \cdot)_K$ 和 $\langle \cdot, \cdot \rangle_e$ 分别为 $K$ 和 $e$ 上的 $L^2$ 内积. 记 $\| \cdot \|_K$ 和 $\| \cdot \|_e$ 分别为 $K$ 和 $e$ 上的 $L^2$ 范数. 简记 $(\cdot, \cdot)$ 和 $\langle \cdot, \cdot \rangle$ 分别为 $\Omega$ 和 $\Gamma$ 上的 $L^2$ 内积, $\| \cdot \| = \| \cdot \|_{L^2(\Omega)}$ . 对于 $\mathcal{M}_h$ 的子集 $T$ 和 $\mathcal{E}_h$ 的子集 $S$ , 分别记

$$
\begin{array}{l} (\cdot , \cdot) _ {\mathcal {T}} = \sum_ {K \in \mathcal {T}} (\cdot , \cdot) _ {K}, \quad \| \cdot \| _ {\mathcal {T}} ^ {2} = \sum_ {K \in \mathcal {T}} \| \cdot \| _ {K} ^ {2}, \\ \langle \cdot , \cdot \rangle_ {\mathcal {S}} = \sum_ {e \in \mathcal {S}} \langle \cdot , \cdot \rangle_ {e}, \quad \| \cdot \| _ {\mathcal {S}} ^ {2} = \sum_ {e \in \mathcal {S}} \| \cdot \| _ {e} ^ {2}. \\ \end{array}
$$

对任意 $e \in \mathcal{E}_h^I$ , 记 $K_1, K_2 \in \mathcal{M}_h$ 为以 $e$ 为公共边的两个单元. 对分片 $H^1$ 函数 $v$ , 记 $v_i = v|_{\partial K_i}$ , 定义 $v$ 在边 $e$ 上的跳量及平均:

$$
[ v ] := v _ {1} - v _ {2}, \quad \{v \} := \frac {v _ {1}}{2} + \frac {v _ {2}}{2}. \tag {8.2}
$$

定义 $\pmb{n}|_{e}$ 为垂直于边 $e$ 的 $\partial K_{1}$ 的单位外法向量.对 $e\in \mathcal{E}_h^D$ ，记 $[v] = \{v\} \coloneqq v,\pmb {n}|_e$ 为垂直于边 $e$ 的单元 $\partial \Omega$ 的单位外法向量.易知，对 $e\in \mathcal{E}_h^I$ 有

$$
[ v w ] = [ v ] \{w \} + \{v \} [ w ]. \tag {8.3}
$$

定义空间

$$
V := \left\{v: v | _ {K} \in H ^ {2} (K), \forall K \in \mathcal {M} _ {h} \right\}. \tag {8.4}
$$

任取 $v \in V$ 乘 (8.1) 的两端, 在 $\Omega$ 上积分, 并在每个单元上利用分部积分公式得

$$
\begin{array}{l} (f, v) = - \int_ {\Omega} \nabla \cdot (\kappa \nabla u) v \mathrm {d} x = - \sum_ {K \in \mathcal {M} _ {h}} \int_ {K} \nabla \cdot (\kappa \nabla u) v \mathrm {d} x \\ = \sum_ {K \in \mathcal {M} _ {h}} \int_ {K} \kappa \nabla u \cdot \nabla v d x - \sum_ {K \in \mathcal {M} _ {h}} \int_ {\partial K} \kappa \nabla u \cdot n _ {K} v d s. \\ \end{array}
$$

由 (8.3) 得

$$
\begin{array}{l} \sum_ {K \in \mathcal {M} _ {h}} \int_ {\partial K} \boldsymbol {\kappa} \nabla u \cdot \boldsymbol {n} _ {K} v \mathrm {d} s = \sum_ {e \in \mathcal {E} _ {h} ^ {I}} \int_ {e} [ \boldsymbol {\kappa} \nabla u \cdot \boldsymbol {n} v ] \mathrm {d} s + \sum_ {e \in \mathcal {E} _ {h} ^ {D}} \int_ {e} \boldsymbol {\kappa} \nabla u \cdot \boldsymbol {n} v \mathrm {d} s \\ = \sum_ {e \in \mathcal {E} _ {h} ^ {I}} \int_ {e} \left(\left[ \boldsymbol {\kappa} \nabla u \cdot \boldsymbol {n} \right] \{v \} + \left\{\boldsymbol {\kappa} \nabla u \cdot \boldsymbol {n} \right\} [ v ]\right) d s + \sum_ {e \in \mathcal {E} _ {h} ^ {D}} \int_ {e} \boldsymbol {\kappa} \nabla u \cdot \boldsymbol {n} v d s \\ = \langle \left\{\kappa \nabla u \cdot \boldsymbol {n} \right\}, [ v ] \rangle_ {\mathcal {E} _ {h} ^ {I D}}, \\ \end{array}
$$

从而精确解 $u$ 满足

$$
\left(\kappa \nabla u, \nabla v\right) _ {\mathcal {M} _ {h}} - \left\langle \left\{\kappa \nabla u \cdot n \right\}, [ v ] \right\rangle_ {\mathcal {E} _ {h} ^ {I D}} = (f, v), \quad \forall v \in V.
$$

定义双线性形式

$$
\begin{array}{l} a _ {h} (u, v) := (\kappa \nabla u, \nabla v) _ {\mathcal {M} _ {h}} - \left(\langle \{\kappa \nabla u \cdot \boldsymbol {n} \}, [ v ] \rangle_ {\mathcal {E} _ {h} ^ {I D}} + \beta \langle [ u ], \{\kappa \nabla v \cdot \boldsymbol {n} \} \rangle_ {\mathcal {E} _ {h} ^ {I D}}\right) + \\ J _ {0} (u, v) + J _ {1} (u, v), \tag {8.5} \\ \end{array}
$$

$$
J _ {0} (u, v) := \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {\gamma_ {0}}{h _ {e}} \int_ {e} [ u ] [ v ] \mathrm {d} s, \tag {8.6}
$$

$$
J _ {1} (u, v) := \sum_ {e \in \mathcal {E} _ {h} ^ {I}} \gamma_ {1} h _ {e} \int_ {e} [ \boldsymbol {\kappa} \nabla u \cdot \boldsymbol {n} ] [ \boldsymbol {\kappa} \nabla v \cdot \boldsymbol {n} ] \mathrm {d} s, \tag {8.7}
$$

其中 $\beta$ 为实数, $\gamma_0 > 0, \gamma_1 \geqslant 0$ 为加罚参数

注8.1 (a) $J_0, J_1$ 称为加罚项，加罚参数 $\gamma_0, \gamma_1$ 可以取得和 $e$ 有关.

(b) $\beta$ 一般取为 $1, -1$ , 或 $0$ . 若 $\beta = 1$ , 则双线性形式 $a_h$ 对称, 否则不对称.

显然对 (8.1) 的精确解 $u$ ，有

$$
\langle [ u ], \{\kappa \nabla v \cdot \boldsymbol {n} \} \rangle_ {\mathcal {E} _ {h} ^ {I D}} = J _ {0} (u, v) = J _ {1} (u, v) = 0, \quad \forall v \in V.
$$

所以（8.1）的精确解 $u$ 满足：

$$
a _ {h} (u, v) = (f, v), \quad \forall v \in V. \tag {8.8}
$$

设 $V_{h}$ 是 $\mathcal{M}_h$ 上的分片线性有限元空间，即

$$
V _ {h} := \left\{v _ {h}: v _ {h} | _ {K} \in P _ {1} (K), \forall K \in \mathcal {M} _ {h} \right\}.
$$

则 $V_{h}\subset V$ .求解(8.1)的内罚间断Galerkin(IPDG)法为：求 $u_{h}\in V_{h}$ 使得

$$
a _ {h} \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in V _ {h}. \tag {8.9}
$$

注8.2 (a) $V_{h}$ 可以换为分片（不连续）高次元空间，甚至可以每个单元用不同次数的多项式.

(b) IPDG 法在 $\beta = 1$ 时简记为 SIPG (Symmetric Interior Penalty Galerkin) 法, 在 $\beta = -1$ 时简记为 NIPG (Nonsymmetric Interior Penalty Galerkin) 法, 在 $\beta = 0$ 时简记为 IIPG (Incomplete Interior Penalty Galerkin) 法.  
(c) 同有限元法比较: IPDG 法更灵活; 容易处理 Dirichlet 边界条件; 满足质量守恒; 但计算量一般较大.

为了进一步的误差分析，引入 $V$ 上的半范数

$$
\| v \| := \left(\left\| \kappa^ {\frac {1}{2}} \nabla v \right\| _ {\mathcal {M} _ {h}} ^ {2} + J _ {0} (v, v) + J _ {1} (v, v) + \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {h _ {e}}{\gamma_ {0}} \| \{\kappa \nabla v \cdot \boldsymbol {n} \} \| _ {L ^ {2} (e)} ^ {2}\right) ^ {\frac {1}{2}}. \tag {8.10}
$$

假设 $\Omega$ 是凸的并且 $\kappa \in C^{1}(\bar{\Omega})^{d\times d}$ 易知 $u\in H^{2}(\varOmega)$

# 8.1.2 对称内罚间断 Galerkin (SIPG) 法的误差分析

本节讨论SIPG法 $(\beta = 1)$ 的 $H^{1}$ 和 $L^2$ 误差估计.先给出双线性形式 $a_{h}$ 的连续性和强制性.

引理8.1

$$
\left| a _ {h} (v, w) \right| \leqslant 2 \| | v \| \| | | w | \|, \quad \forall v, w \in V. \tag {8.11}
$$

存在常数 $\underline{\gamma} > 0$ 与 $h$ 及加罚参数无关, 使得当 $\gamma_0 \geqslant \underline{\gamma}$ 时,

$$
a _ {h} \left(v _ {h}, v _ {h}\right) \geqslant \frac {1}{2} \| \left| v _ {h} \right| \| ^ {2}, \quad \forall v _ {h} \in V _ {h}. \tag {8.12}
$$

证明 由Cauchy-Schwarz不等式易知(8.11)成立.下面证(8.12).由(8.5)—(8.7)和(8.10),

$$
\begin{array}{l} a _ {h} (v _ {h}, v _ {h}) = \| \left| v _ {h} \right| \| ^ {2} - 2 \langle \{\kappa \nabla v _ {h} \cdot \boldsymbol {n} \}, [ v _ {h} ] \rangle_ {\mathcal {E} _ {h} ^ {I D}} - \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {h _ {e}}{\gamma_ {0}} \| \{\kappa \nabla v _ {h} \cdot \boldsymbol {n} \} \| _ {e} ^ {2} \\ \geqslant \| v _ {h} \| ^ {2} - \frac {1}{2} \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {\gamma_ {0}}{h _ {e}} \int_ {e} [ v _ {h} ] ^ {2} - 3 \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {h _ {e}}{\gamma_ {0}} \| \left\{\kappa \nabla v _ {h} \cdot n \right\} \| _ {e} ^ {2}. \\ \end{array}
$$

由局部的迹不等式 (6.77) 得

$$
\begin{array}{l} 3 \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {h _ {e}}{\gamma_ {0}} \| \left\{\kappa \nabla v _ {h} \cdot \boldsymbol {n} \right\} \| _ {e} ^ {2} \leqslant C \sum_ {K \in \mathcal {M} _ {h}} \frac {h _ {K}}{\gamma_ {0}} h _ {K} ^ {- 1} \| \nabla v _ {h} \| _ {K} ^ {2} \\ \leqslant \frac {C}{\gamma_ {0}} \sum_ {K \in \mathcal {M} _ {h}} \left\| \kappa^ {\frac {1}{2}} \nabla v _ {h} \right\| _ {K} ^ {2}. \\ \end{array}
$$

从而

$$
a _ {h} \left(v _ {h}, v _ {h}\right) \geqslant \| v _ {h} \| ^ {2} - \max  \left\{\frac {1}{2}, \frac {C}{\gamma_ {0}} \right\} \| v _ {h} \| ^ {2}.
$$

所以 $\gamma_0$ 足够大时 (8.12) 成立. 证毕

然后给出 Céa 引理

引理8.2 在引理8.1的条件下，

$$
\| u - u _ {h} \| \lesssim \inf  _ {v _ {h} \in V _ {h}} \| u - v _ {h} \|,
$$

其中，“ $A \lesssim B$ ”表示存在不依赖步长 $h$ 的正常数 $C$ ，使得 $A \leqslant CB$ .

证明 由 (8.8) 和 (8.9) 得：

$$
a _ {h} \left(u - u _ {h}, v _ {h}\right) = 0, \quad \forall v _ {h} \in V _ {h}. \tag {8.13}
$$

故由引理8.1,

$$
\begin{array}{l} \left\| \left| u _ {h} - v _ {h} \right| \right\| ^ {2} \leqslant 2 a _ {h} (u _ {h} - v _ {h}, u _ {h} - v _ {h}) = 2 a _ {h} (u - v _ {h}, u _ {h} - v _ {h}) \\ \leqslant 4 \| u - u _ {h} \| \| u _ {h} - v _ {h} \|. \\ \end{array}
$$

得

$$
\| \left. u _ {h} - v _ {h} \right\| \leqslant 4 \| \left. u - v _ {h} \right\|.
$$

再由三角不等式知

$$
\| u - u _ {h} \| \leqslant \| u - v _ {h} \| + \| v _ {h} - u _ {h} \| \leqslant 5 \| u - v _ {h} \|.
$$

证毕.

再给出插值估计. 记 $I_h u$ 为 $u$ 的有限元插值

引理8.3 设 $u \in H^2(\Omega)$ , 则

$$
\| u - I _ {h} u \| \lesssim h ^ {2} | u | _ {H ^ {2} (\Omega)}, \quad \| | u - I _ {h} u | \| \lesssim h \left(1 + \gamma_ {0} + \gamma_ {0} ^ {- 1} + \gamma_ {1}\right) ^ {\frac {1}{2}} | u | _ {H ^ {2} (\Omega)}.
$$

证明 记 $I_{K}u$ 为 $\mathcal{U}$ 在单元 $K$ 上的有限元插值.我们有

$$
\left\| u - I _ {K} u \right\| _ {L ^ {2} (K)} + h _ {K} \| u - I _ {K} u \| _ {H ^ {1} (K)} \lesssim h _ {K} ^ {2} | u | _ {H ^ {2} (K)}.
$$

则 $I_{h}u$ 的第一个估计显然成立.下面证第二个.记 $\eta_h = u - I_hu.$ 由（8.10）及局部迹不等式，

$$
\begin{array}{l} \| \eta_ {h} \| ^ {2} = \left\| \kappa^ {\frac {1}{2}} \nabla \eta_ {h} \right\| ^ {2} + \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {\gamma_ {0}}{h _ {e}} \| [ \eta_ {h} ] \| _ {e} ^ {2} + \\ \sum_ {e \in \mathcal {E} _ {h} ^ {I}} \gamma_ {1} h _ {e} \| [ \kappa \nabla \eta_ {h} \cdot \boldsymbol {n} ] \| _ {e} ^ {2} + \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {h _ {e}}{\gamma_ {0}} \| \{\kappa \nabla \eta_ {h} \cdot \boldsymbol {n} \} \| _ {e} ^ {2} \\ \lesssim \| \nabla \eta_ {h} \| ^ {2} + \sum_ {K \in \mathcal {M} _ {h}} \frac {\gamma_ {0}}{h} \left(h ^ {- 1} \| \eta_ {h} \| _ {K} ^ {2} + h \| \nabla \eta_ {h} \| _ {K} ^ {2}\right) + \\ \sum_ {K \in \mathcal {M} _ {h}} \left(\gamma_ {1} + \gamma_ {0} ^ {- 1}\right) h \left(h ^ {- 1} \| \nabla \eta_ {h} \| _ {K} ^ {2} + h | \nabla \eta_ {h} | _ {H ^ {1} (K)} ^ {2}\right) \\ \lesssim \left(1 + \gamma_ {0} + \gamma_ {0} ^ {- 1} + \gamma_ {1}\right) \left(h ^ {- 2} \| \eta_ {h} \| ^ {2} + \| \nabla \eta_ {h} \| ^ {2} + h ^ {2} | \nabla u | _ {H ^ {1} (\Omega)} ^ {2}\right) \\ \lesssim h ^ {2} \left(1 + \gamma_ {0} + \gamma_ {0} ^ {- 1} + \gamma_ {1}\right) | u | _ {H ^ {2} (\Omega)} ^ {2}. \\ \end{array}
$$

证毕.

注8.3 若网格 $\mathcal{M}_h$ 是协调的，则 $I_h u \in V_h \cap H^1(\Omega)$ ，且在 $\Gamma_D$ 上 $I_h u = 0$ 又 $\|[\eta_h]\|_e^2 = 0$ ，从而第二个估计可改进为

$$
\| \| u - I _ {h} u \| \lesssim \left(1 + \gamma_ {0} ^ {- 1} + \gamma_ {1}\right) ^ {\frac {1}{2}} | u | _ {H ^ {2} (\Omega)}.
$$

结合引理8.2和8.3可得 $H^{1}$ 误差估计

定理8.1 在引理8.1的条件下，

$$
\| \| u - u _ {h} \| \| \lesssim h (1 + \gamma_ {0} + \gamma_ {1}) ^ {\frac {1}{2}} | u | _ {H ^ {2} (\Omega)}.
$$

注8.4 若网格 $\mathcal{M}_h$ 是协调的，则上面估计可改进为

$$
\left\| u - u _ {h} \right\| \leqslant C h \left(1 + \gamma_ {1}\right) ^ {\frac {1}{2}} | u | _ {H ^ {2} (\Omega)}.
$$

下面讨论 $L^2$ 误差估计. 考虑对偶问题:

$$
- \nabla \cdot (\kappa \nabla w) = u - u _ {h}, \quad \boldsymbol {x} \in \Omega , \quad w | _ {\Gamma_ {D}} = 0, \quad a \nabla w \cdot \boldsymbol {n} | _ {\Gamma_ {N}} = 0. \tag {8.14}
$$

易知 $w\in H^{2}(\varOmega)$ 且

$$
\left\| w \right\| _ {H ^ {2} (\Omega)} \lesssim \left\| u - u _ {h} \right\|.
$$

类似于 (8.8) 可知 $w$ 满足

$$
a _ {h} (w, v) = (u - u _ {h}, v), \quad \forall v \in V. \tag {8.15}
$$

取 $v = u - u_{h}$ 并利用(8.13）及引理8.1—8.3得

$$
\begin{array}{l} \| u - u _ {h} \| _ {L ^ {2} (\Omega)} ^ {2} = a _ {h} (w, u - u _ {h}) = a _ {h} (u - u _ {h}, w) \\ = a _ {h} \left(u - u _ {h}, w - I _ {h} w\right) \leqslant C \| \left\| u - u _ {h} \right\| \| \left\| w - I _ {h} \right\| \\ \lesssim \| u - u _ {h} \| h \left(1 + \gamma_ {0} + \gamma_ {0} ^ {- 1} + \gamma_ {1}\right) ^ {\frac {1}{2}} | w | _ {H ^ {2} (\Omega)} \\ \lesssim h \left(1 + \gamma_ {0} + \gamma_ {0} ^ {- 1} + \gamma_ {1}\right) ^ {\frac {1}{2}} \| u - u _ {h} \| \| u - u _ {h} \|. \\ \end{array}
$$

最后由定理8.1可得如下 $L^2$ 误差估计：

定理8.2 在引理8.1的条件下，

$$
\| u - u _ {h} \| \lesssim h ^ {2} (1 + \gamma_ {0} + \gamma_ {1}) | u | _ {H ^ {2} (\Omega)}.
$$

注8.5 若网格 $\mathcal{M}_h$ 是协调的，则 $L^2$ 估计可改进为

$$
\left\| u - u _ {h} \right\| \leqslant C h ^ {2} (1 + \gamma_ {1}) | u | _ {H ^ {2} (\Omega)}.
$$

# 8.1.3 非对称内罚间断 Galerkin (NIPG) 法的误差分析

本节讨论NIPG法 $(\beta = -1)$ 的 $H^{1}$ 和 $L^2$ 误差估计

先给出双线性形式 $a_{h}$ 的连续性和强制性

引理8.4

$$
\left| a _ {h} (v, w) \right|, \left| a _ {h} (w, v) \right| \leqslant 2 \| | v \| \| | | w | \|, \quad \forall v, w \in V. \tag {8.16}
$$

$$
a _ {h} \left(v _ {h}, v _ {h}\right) \geqslant \frac {\gamma_ {0}}{\gamma_ {0} + \alpha} \| v _ {h} \| ^ {2}, \quad \forall v _ {h} \in V _ {h}, \tag {8.17}
$$

其中常数 $\alpha > 0$ 与 $h$ 及加罚参数无关.

证明 由Cauchy-Schwarz不等式易知(8.16)成立.下面证(8.17).显然

$$
a _ {h} (v _ {h}, v _ {h}) = \| v _ {h} \| ^ {2} - \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {h _ {e}}{\gamma_ {0}} \| \kappa \nabla v _ {h} \cdot \boldsymbol {n} \| _ {e} ^ {2}, \quad \forall v _ {h} \in V _ {h}.
$$

另外，由引理8.1的证明知

$$
\cdot \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {h _ {e}}{\gamma_ {0}} \| \left\{\kappa \nabla v _ {h} \cdot \boldsymbol {n} \right\} \| _ {e} ^ {2} \leqslant \frac {C}{\gamma_ {0}} \sum_ {K \in \mathcal {M} _ {h}} \left\| \kappa^ {\frac {1}{2}} \nabla v _ {h} \right\| _ {K} ^ {2} \leqslant \frac {C}{\gamma_ {0}} a _ {h} (v _ {h}, v _ {h}).
$$

故

$$
\| v _ {h} \| ^ {2} \leqslant \left(1 + \frac {C}{\gamma_ {0}}\right) a _ {h} (v _ {h}, v _ {h}).
$$

证毕.

类似于引理8.2可得NIPG法的Céa引理，证明略去

引理8.5 假设 $\gamma_0\gtrsim 1$ ，则

$$
\| u - u _ {h} \| \lesssim \inf  _ {v _ {h} \in V _ {h}} \| u - v _ {h} \|.
$$

其中，“ $A\gtrsim B$ ”表示存在不依赖步长 $h$ 的正常数 $C$ ，使得 $A\geqslant CB$

结合引理8.3及8.5可得NIPG法的 $H^{1}$ 误差估计

定理8.3 假设 $\gamma_0\gtrsim 1$ ，则

$$
\| u - u _ {h} \| \lesssim h (1 + \gamma_ {0} + \gamma_ {1}) ^ {\frac {1}{2}} | u | _ {H ^ {2} (\Omega)}.
$$

注8.6 (a) 若网格 $\mathcal{M}_h$ 是协调的, 则上面估计可改进为

$$
\left\| u - u _ {h} \right\| \lesssim h \left(1 + \gamma_ {1}\right) ^ {\frac {1}{2}} | u | _ {H ^ {2} (\Omega)}.
$$

(b) 同SIPG法比，NIPG法不要求 $\gamma_0$ 足够大，但刚度矩阵不对称.

下面讨论 $L^2$ 误差估计. 考虑对偶问题. 我们仍然有

$$
\| w \| _ {H ^ {2} (\Omega)} \lesssim \| u - u _ {h} \| _ {L ^ {2} (\Omega)}.
$$

$$
a _ {h} (w, v) = (u - u _ {h}, v), \quad \forall v \in V.
$$

取 $v = u - u_{h}$ ，注意到 $a_{h}$ 不对称，并利用(8.13）及引理8.1—8.3得

$$
\begin{array}{l} \| u - u _ {h} \| _ {L ^ {2} (\Omega)} ^ {2} = a _ {h} (w, u - u _ {h}) \\ = a _ {h} (u - u _ {h}, w) - 2 \langle [ u - u _ {h} ], \kappa \nabla w \cdot n \rangle_ {\mathcal {E} _ {h} ^ {I D}} \\ = a _ {h} \left(u - u _ {h}, w - I _ {h} w\right) - 2 \langle [ u - u _ {h} ], \kappa \nabla w \cdot n \rangle_ {\mathcal {E} _ {h} ^ {I D}} \\ \lesssim \| u - u _ {h} \| \| w - I _ {h} w \| + \\ J _ {0} (u - u _ {h}, u - u _ {h}) ^ {\frac {1}{2}} \left(\sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {h _ {e}}{\gamma_ {0}} \| \kappa \nabla w \cdot n \| _ {e} ^ {2}\right) ^ {\frac {1}{2}} \\ \lesssim \| u - u _ {h} \| \left(h \left(1 + \gamma_ {0} + \gamma_ {0} ^ {- 1} + \gamma_ {1}\right) ^ {\frac {1}{2}} | w | _ {H ^ {2} (\Omega)} + \right. \\ \gamma_ {0} ^ {- \frac {1}{2}} \left(\left\| \nabla w \right\| + h | \nabla w | _ {H ^ {1} (\Omega)}\right) \Bigg) \\ \lesssim \left(h \left(1 + \gamma_ {0} + \gamma_ {0} ^ {- 1} + \gamma_ {1}\right) ^ {\frac {1}{2}} + \gamma_ {0} ^ {- \frac {1}{2}}\right) \| | u - u _ {h} | \| | | u - u _ {h} |. \\ \end{array}
$$

最后由定理8.3可得如下 $L^2$ 误差估计：

定理8.4 假设 $\gamma_0\gtrsim 1$ ，则

$$
\left\| u - u _ {h} \right\| \lesssim \left(h ^ {2} + \gamma_ {0} ^ {- \frac {1}{2}} h\right) \left(1 + \gamma_ {0} + \gamma_ {1}\right) | u | _ {H ^ {2} (\Omega)}.
$$

注8.7 (a)上面NIPG方法的 $L^2$ 估计只是一阶收敛

(b) 若网格 $\mathcal{M}_h$ 是协调的, 则 $L^2$ 估计可改进为

$$
\left\| u - u _ {h} \right\| \lesssim \left(h ^ {2} + \gamma_ {0} ^ {- \frac {1}{2}} h\right) (1 + \gamma_ {1}) | u | _ {H ^ {2} (\Omega)}.
$$

所以，若 $\gamma_0\gtrsim h^{-2}$ ，则

$$
\| u - u _ {h} \| \lesssim h ^ {2} (1 + \gamma_ {1}) | u | _ {H ^ {2} (\Omega)}.
$$

但此时刚度矩阵的条件数较大，为 $O(h^{-4})$

(c) 对高次元的 NIPG 方法, 通过仔细估计 $\|\nabla w\|$ , 可以证明 $L^2$ 误差可以达到满阶收敛.

# 8.2 局部间断 Galerkin (LDG) 法

LDG法[38]在计算流体等领域有重要应用

# 8.2.1 离散格式

首先引入中间变量 $\sigma = \kappa \delta$ ， $\delta = \nabla u$ 将椭圆问题(8.1)改写为如下一阶偏微分方程组：

$$
\sigma = \kappa \delta , \delta = \nabla u, - \nabla \cdot \sigma = f, \quad \boldsymbol {x} \in \Omega , \quad u | _ {\Gamma_ {D}} = 0, \quad (\sigma \cdot n) | _ {\Gamma_ {N}} = 0. \tag {8.18}
$$

用检验函数 $\theta, \tau$ 和 $v$ 分别乘前三个方程，并在任一单元 $K \in \mathcal{M}_h$ 上积分得：

$$
\begin{array}{l} \int_ {K} \sigma \cdot \theta d x = \int_ {K} \kappa \delta \cdot \theta d x. \\ \int_ {K} \delta \cdot \tau \mathrm {d} \boldsymbol {x} = - \int_ {K} u \nabla \cdot \tau \mathrm {d} \boldsymbol {x} + \int_ {\partial K} u \tau \cdot n _ {K} \mathrm {d} s, \\ \int_ {K} \sigma \cdot \nabla v \mathrm {d} x = \int_ {K} f v \mathrm {d} x + \int_ {\partial K} \sigma \cdot n _ {K} v \mathrm {d} s. \\ \end{array}
$$

为了离散化上面的变分形式，定义如下的逼近空间：

$$
V _ {h} := \left\{v _ {h} \in L ^ {2} (\Omega): v _ {h} | _ {K} \in P _ {1} (K), \forall K \in \mathcal {M} _ {h} \right\}, \quad \Sigma_ {h} := \left(V _ {h}\right) ^ {d}. \tag {8.19}
$$

则求解椭圆型问题 (8.1) 的LDG方法为: 求 $u_{h} \in V_{h}, \sigma_{h}, \delta_{h} \in \Sigma_{h}$ 使得

$$
\int_ {K} \sigma_ {h} \cdot \theta_ {h} \mathrm {d} \boldsymbol {x} = \int_ {K} \kappa \delta_ {h} \cdot \theta_ {h} \mathrm {d} \boldsymbol {x}, \quad \forall \theta_ {h} \in \Sigma_ {h}, \tag {8.20}
$$

$$
\int_ {K} \delta_ {h} \cdot \tau_ {h} \mathrm {d} \boldsymbol {x} = - \int_ {K} u _ {h} \nabla \cdot \tau_ {h} \mathrm {d} \boldsymbol {x} + \int_ {\partial K} \hat {u} _ {h} \tau_ {h} \cdot n _ {K} \mathrm {d} s, \quad \forall \tau_ {h} \in \Sigma_ {h}, \tag {8.21}
$$

$$
\int_ {K} \sigma_ {h} \cdot \nabla v _ {h} d x = \int_ {K} f v _ {h} d x + \int_ {\partial K} \hat {\sigma} _ {h} \cdot n _ {K} v _ {h} d s, \quad \forall v _ {h} \in V _ {h}, K \in \mathcal {M} _ {h}. \tag {8.22}
$$

其中 $\hat{\sigma}_h$ 和 $\hat{u}_h$ 称为数值流通量，分别是 $\kappa \nabla u$ 和 $u$ 的近似.为了定义数值流通量，我们引入如下跳量的定义.对内部边 $e = K_{1}\cap K_{2}\in \mathcal{E}_{h}^{I}$ ，定义

$$
\llbracket \varphi \rrbracket := \varphi_ {1} \cdot n _ {K _ {1}} + \varphi_ {2} \cdot n _ {K _ {2}}.
$$

显然, 如果 $\varphi$ 是标量函数, 则 $[v]$ 是向量并与边 $e$ 垂直; 如果 $\varphi$ 是向量函数, 则 $[v]$ 是标量. 另外, 由 (8.2) 知, 两种跳量关系为 $[\varphi] = [\varphi \cdot n_{K_1}]$ . 对应 $\Gamma$ 上的边 $e$ , 规定

$$
\llbracket \varphi \rrbracket := \varphi \cdot n.
$$

其中 $n$ 是 $\partial \Omega$ 单位外法向量. 下面给出 (8.20)一(8.22) 中数值流通量在边 $e$ 上的定义:

$$
\hat {u} _ {h} \mid_ {e}: = \left\{ \begin{array}{l l} \left\{u _ {h} \right\} - \beta_ {e} \cdot \llbracket u _ {h} \rrbracket , & e \in \mathcal {E} _ {h} ^ {I}, \\ 0, & e \in \mathcal {E} _ {h} ^ {D}, \\ u _ {h}, & e \in \mathcal {E} _ {h} ^ {N}. \end{array} \right. \tag {8.23}
$$

$$
\hat {\sigma} _ {h} \mid_ {e}: = \left\{ \begin{array}{l l} \left\{\sigma_ {h} \right\} + \beta_ {e} \llbracket \sigma_ {h} \rrbracket - \frac {\gamma_ {e}}{h _ {e}} \llbracket u _ {h} \rrbracket , & e \in \mathcal {E} _ {h} ^ {I}, \\ \sigma_ {h} - \frac {\gamma_ {e}}{h _ {e}} \llbracket u _ {h} \rrbracket , & e \in \mathcal {E} _ {h} ^ {D}, \\ 0, & e \in \mathcal {E} _ {h} ^ {N}. \end{array} \right. \tag {8.24}
$$

其中 $\beta_{e}\in \mathbb{R}^{d},\gamma_{e}\in \mathbb{R}^{+}$ .LDG格式由(8.20)—(8.24）组成.如果系数 $\kappa$ 是分片常数的，则显然 $\sigma_h = \kappa \delta_h$ ，此时可以不引入中间变量 $\delta_h$

# 8.2.2 原始变量形式

为了进行误差估计，我们将LDG法(8.20)一(8.24）中的 $\delta_h,\sigma_h$ 消掉，改写为关于原始变量 $u_{h}$ 的公式.记 $\nabla_{h}$ 为 $\mathcal{M}_h$ 上的分片梯度算子，即

$$
\left. \left(\nabla_ {h} v\right) \right| _ {K} = \left. \nabla \left(v \right| _ {K}\right), \quad \left. \left(\nabla_ {h} \cdot \tau\right) \right| _ {K} = \nabla \cdot \left(\tau \right| _ {K}), \quad \forall K \in \mathcal {M} _ {h}.
$$

将 (8.21)一(8.22) 按 $K \in \mathcal{M}_h$ 求和得

$$
\begin{array}{l} \int_ {\Omega} \delta_ {h} \cdot \tau_ {h} d x = - \int_ {\Omega} u _ {h} \nabla_ {h} \cdot \tau_ {h} d x + \sum_ {K \in \mathcal {M} _ {h}} \int_ {\partial K} \hat {u} _ {h} \tau_ {h} \cdot n _ {K} d x \\ = \int_ {\Omega} \nabla_ {h} u _ {h} \cdot \tau_ {h} \mathrm {d} \boldsymbol {x} + \sum_ {K \in \mathcal {M} _ {h}} \int_ {\partial K} \left(\hat {u} _ {h} - u _ {h}\right) \tau_ {h} \cdot n _ {K} \mathrm {d} \boldsymbol {x}, \\ \end{array}
$$

$$
\int_ {\Omega} \sigma_ {h} \cdot \nabla_ {h} v _ {h} = \int_ {\Omega} f v _ {h} d x + \sum_ {K \in \mathcal {M} _ {h}} \int_ {\partial K} \hat {\sigma} _ {h} \cdot n _ {K} v _ {h} d x.
$$

对于标量函数 $v \in \Pi_{K \in \mathcal{M}_h} L^2(\partial K)$ 和向量值函数 $\varphi \in \Pi_{K \in \mathcal{M}_h} (L^2(\partial K))^d$ , 考虑和式 $\sum_{K \in \mathcal{M}_h} \int_{\partial K} v \varphi \cdot n_K \mathrm{d}s.$ 易知

$$
\begin{array}{l} \sum_ {K \in \mathcal {M} _ {h}} \int_ {\partial K} v \varphi \cdot n _ {K} d s = \sum_ {e \in \mathcal {E} _ {h} ^ {D} \cup \mathcal {E} _ {h} ^ {N}} \int_ {e} v \varphi \cdot n d s + \sum_ {e \in \mathcal {E} _ {h} ^ {I}} \int_ {e} [   [ v \varphi ]   ] d s \\ = \langle [   [ v ] ], \{\varphi \} \rangle_ {\mathcal {E} _ {h} ^ {I D}} + \langle \{v \}, [   [ \varphi ] ] \rangle_ {\mathcal {E} _ {h} ^ {I N}}. \tag {8.25} \\ \end{array}
$$

记

$$
\beta | _ {e} = \left\{ \begin{array}{l l} \beta_ {e}, & e \in \mathcal {E} _ {h} ^ {I}, \\ 0, & e \in \mathcal {E} _ {h} ^ {D}. \end{array} \right.
$$

由以上三式及 (8.23)一(8.24)得

$$
\begin{array}{l} \int_ {\Omega} \delta_ {h} \cdot \tau_ {h} d x = \int_ {\Omega} \nabla_ {h} u _ {h} \cdot \tau_ {h} d x - \langle [   [ u _ {h} ]   ], \{\tau_ {h} \} \rangle_ {\mathcal {E} _ {h} ^ {I D}} + \langle \hat {u} _ {h} - \{u _ {h} \}, [   [ \tau_ {h} ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I}} \\ = \int_ {\Omega} \nabla_ {h} u _ {h} \cdot \tau_ {h} d x - \langle [   [ u _ {h} ] ], \{\tau_ {h} \} + \beta [   [ \tau_ {h} ] ] \rangle_ {\mathcal {E} _ {h} ^ {I D}}. \tag {8.26} \\ \end{array}
$$

$$
\begin{array}{l} \int_ {\Omega} \sigma_ {h} \cdot \nabla_ {h} v _ {h} d \boldsymbol {x} = \int_ {\Omega} f v _ {h} d \boldsymbol {x} + \langle [   [ v _ {h} ] ], \hat {\sigma} _ {h} \rangle_ {\mathcal {E} _ {h} ^ {I D}} \\ = \int_ {\Omega} f v _ {h} \mathrm {d} \boldsymbol {x} + \left\langle \llbracket v _ {h} \rrbracket , \left\{\sigma_ {h} \right\} + \beta \llbracket \sigma_ {h} \rrbracket \right\rangle_ {\mathcal {E} _ {h} ^ {I D}} - \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {\gamma_ {e}}{h _ {e}} \int_ {e} \llbracket u _ {h} \rrbracket \cdot \llbracket v _ {h} \rrbracket \mathrm {d} s. \tag {8.27} \\ \end{array}
$$

引入提升算子 $L_{h}:V_{h} + H^{1}(\varOmega)\mapsto \Sigma_{h}$

$$
\int_ {\Omega} L _ {h} v \cdot \varphi_ {h} \mathrm {d} x = \langle [   [ v ] ], \{\varphi_ {h} \} + \beta [   [ \varphi_ {h} ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I D}}, \quad \forall \varphi_ {h} \in \Sigma_ {h}. \tag {8.28}
$$

由 (8.26) 得

$$
\delta_ {h} = \nabla_ {h} u _ {h} - L _ {h} u _ {h}. \tag {8.29}
$$

另外，由(8.20)对 $K$ 求和得：

$$
\int_ {\Omega} \sigma_ {h} \cdot \theta_ {h} \mathrm {d} \boldsymbol {x} = \int_ {\Omega} \kappa \delta_ {h} \cdot \theta_ {h} \mathrm {d} \boldsymbol {x}, \quad \forall \theta_ {h} \in \Sigma_ {h}. \tag {8.30}
$$

由 (8.27)一(8.30),

$$
\begin{array}{l} \int_ {\Omega} f v _ {h} \mathrm {d} \boldsymbol {x} = \int_ {\Omega} \sigma_ {h} \cdot \nabla_ {h} v _ {h} \mathrm {d} \boldsymbol {x} - \int_ {\Omega} \sigma_ {h} \cdot L _ {h} v _ {h} \mathrm {d} \boldsymbol {x} + \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {\gamma_ {e}}{h _ {e}} \int_ {e} [   [ u _ {h} ]   ] \cdot [   [ v _ {h} ]   ] \mathrm {d} s \\ = \int_ {\Omega} \kappa \delta_ {h} \cdot (\nabla_ {h} v _ {h} - L _ {h} v _ {h}) d \boldsymbol {x} + \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {\gamma_ {e}}{h _ {e}} \int_ {e} [   [ u _ {h} ]   ] \cdot [   [ v _ {h} ]   ] d s \\ \end{array}
$$

$$
= \int_ {\Omega} \kappa (\nabla_ {h} u _ {h} - L _ {h} u _ {h}) \cdot (\nabla_ {h} v _ {h} - L _ {h} v _ {h}) d \boldsymbol {x} + \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {\gamma_ {e}}{h _ {e}} \int_ {e} [   [ u _ {h} ]   ] \cdot [   [ v _ {h} ]   ] d s.
$$

引入 $V_{h} + H^{1}(\Omega)$ 上的双线性形式

$$
B _ {h} (u, v) := \int_ {\Omega} \kappa (\nabla_ {h} u - L _ {h} u) \cdot (\nabla_ {h} v - L _ {h} v) d x + \sum_ {e \in \mathcal {E} _ {h} ^ {I D}} \frac {\gamma_ {e}}{h _ {e}} \int_ {e} [   [ u ]   ] \cdot [   [ v ]   ] d s. \tag {8.31}
$$

得到LDG的原始变量公式(primalformulation):求 $u_{h}\in V_{h}$ 使得

$$
B _ {h} \left(u _ {h}, v _ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in V _ {h}. \tag {8.32}
$$

# 8.2.3 误差估计

不妨设 $\gamma_{e} \equiv \gamma > 0$ . 定义离散的能量范数

$$
\| v \| _ {h} ^ {2} = \left\| \kappa^ {\frac {1}{2}} \nabla_ {h} v \right\| ^ {2} + \gamma \left\| h ^ {- \frac {1}{2}} [ v ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} ^ {2}, \tag {8.33}
$$

其中 $\left\| h^{-\frac{1}{2}}[[v]]\right\|_{\mathcal{E}_h^{ID}}^2 = \sum_{e\in \mathcal{E}_h^{ID}}h_e^{-1}\| [[v]]\| _e^2.$

$B_{h}$ 的连续性和强制性

先给出 $L_{h}$ 的稳定性估计

引理8.6 设 $|\beta_{e}| \lesssim 1$ ，则存在常数 $C_{L} > 0$ 使得

$$
\left\| \kappa^ {\frac {1}{2}} L _ {h} v \right\| \leqslant C _ {L} \left\| h ^ {- \frac {1}{2}} [ v ] \right\| _ {\mathcal {E} _ {h} ^ {I D}}, \quad \forall v \in V _ {h} + H ^ {1} (\Omega).
$$

证明 由（8.28）及局部迹不等式和逆估计得

$$
\begin{array}{l} \int_ {\Omega} L _ {h} v \cdot \varphi_ {h} \mathrm {d} x \lesssim \left(\sum_ {e \in \mathcal {E} _ {h} ^ {I D}} h _ {e} ^ {- 1} \| [   [ v ]   ] \| _ {e} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {e \in \mathcal {E} _ {h} ^ {I D}} h _ {e} \left(\| \{\varphi_ {h} \} \| _ {e} + \| [   [ \varphi_ {h} ]   ] \| _ {e}\right) ^ {2}\right) ^ {\frac {1}{2}} \\ \lesssim \left(\sum_ {e \in \mathcal {E} _ {h} ^ {I D}} h _ {e} ^ {- 1} \| [ v ] \| _ {e} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {K \in \mathcal {M} _ {h}} h _ {K} \| \varphi_ {h} \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \\ \lesssim \left\| h ^ {- \frac {1}{2}} [   [ v ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} \| \varphi_ {h} \|, \\ \end{array}
$$

取 $\varphi_h = L_h v$ 即得证明

下面引理给出双线性形式 $B_{h}$ 的连续性和强制性

引理8.7 设 $|\beta_{e}| \lesssim 1$ ，则

$$
B _ {h} (u, v) \lesssim \frac {1 + \gamma}{\gamma} \| u \| _ {h} \| v \| _ {h}, \quad \forall u, v \in V _ {h} + H ^ {1} (\Omega), \tag {8.34}
$$

$$
B _ {h} (v, v) \gtrsim \frac {\gamma}{1 + \gamma} \| v \| _ {h} ^ {2}, \quad \forall v \in V _ {h} + H ^ {1} (\Omega). \tag {8.35}
$$

证明 由 (8.31) (8.33) 及引理 8.6 得

$$
\begin{array}{l} B _ {h} (u, v) \leqslant \left(\| \kappa^ {\frac {1}{2}} \nabla_ {h} u \| + \| \kappa^ {\frac {1}{2}} L _ {h} u \|\right) \left(\| \kappa^ {\frac {1}{2}} \nabla_ {h} v \| + \| \kappa^ {\frac {1}{2}} L _ {h} v \|\right) + \\ \gamma \left\| h ^ {- \frac {1}{2}} [   [ u ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} \left\| h ^ {- \frac {1}{2}} [   [ v ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} \\ \leqslant \left[ \left(1 + C _ {L} \gamma^ {- \frac {1}{2}}\right) ^ {2} + 1 \right] \| | u | | | h | | | v | | | _ {h}, \\ \end{array}
$$

即 (8.34) 成立. 下面证明 (8.35). 对任意 $\varepsilon \in (0,1)$ ,

$$
\begin{array}{l} B _ {h} (v, v) = \left\| \kappa^ {\frac {1}{2}} \left(\nabla_ {h} v - L _ {h} v\right) \right\| _ {L ^ {2} (\Omega)} ^ {2} + \gamma \left\| h ^ {- \frac {1}{2}} [   [ v ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} ^ {2} \\ \geqslant (1 - \varepsilon) \| \kappa^ {\frac {1}{2}} \nabla_ {h} v \| ^ {2} + (1 - \varepsilon^ {- 1}) \| \kappa^ {\frac {1}{2}} L _ {h} v \| ^ {2} + \gamma \| h ^ {- \frac {1}{2}} [   [ v ]   ] \| _ {\varepsilon_ {h} ^ {I D}} ^ {2} \\ \geqslant (1 - \varepsilon) \left\| \kappa^ {\frac {1}{2}} \nabla_ {h} v \right\| ^ {2} + (\gamma - C _ {L} ^ {2} (\varepsilon^ {- 1} - 1)) \left\| h ^ {- \frac {1}{2}} [   [ v ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} ^ {2}. \\ \end{array}
$$

取 $\varepsilon = \frac{C_L^2}{C_L^2 + \gamma / 2}$ 即知（8.35）成立.证毕

相容性

下面考虑 Galerkin 正交性. 设 $u$ 和 $u_h$ 分别是椭圆型问题 (8.1) 及其 LDG 离散 (8.20) — (8.24) 的解. 记 $Q_h$ 为 $L^2(\Omega)^d$ 到 $\Sigma_h$ 的正交 $L^2$ 投影. 由 (8.28) 易知 $L_h u = 0$ , 因此, 由 (8.31) 得: $\forall v \in V_h + H^1(\Omega)$ ,

$$
\begin{array}{l} B _ {h} (u, v) = \int_ {\Omega} \kappa \nabla u \cdot (\nabla_ {h} v - L _ {h} v) d x \\ = \int_ {\Omega} \kappa \nabla u \cdot \nabla_ {h} v \mathrm {d} x - \int_ {\Omega} Q _ {h} (\kappa \nabla u) \cdot L _ {h} v \mathrm {d} x \\ = - \int_ {\Omega} \nabla \cdot (\kappa \nabla u) v \mathrm {d} x + \sum_ {K \in \mathcal {M} _ {h}} \int_ {\partial K} \kappa \nabla u \cdot n _ {K} v \mathrm {d} s - \\ \langle [   [ v ] ], \{Q _ {h} (\kappa \nabla u) \} + \beta [   [ Q _ {h} (\kappa \nabla u) ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I D}} \\ = (f, v) + \langle [   [ v ]   ], \{\kappa \nabla u - Q _ {h} (\kappa \nabla u) \} + \beta [   [ \kappa \nabla u - Q _ {h} (\kappa \nabla u) ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I D}}. \\ \end{array}
$$

可以看出， $u_{h}$ 不是精确满足Galerkin正交性，即 $B_{h}(u - u_{h},v_{h})$ 可能不为零

定义残量

$$
R _ {h} (u, v) := B _ {h} (u, v) - (f, v). \tag {8.36}
$$

显然， $B_{h}(u - u_{h},v_{h}) = R_{h}(u,v_{h})$ .下面引理给出残量 $R_{h}(u,v)$ 的估计

引理8.8 设 $|\beta_{e}|\lesssim 1,u$ 是椭圆型问题(8.1)的解，则

$$
\left| R _ {h} (u, v) \right| \lesssim h \| \nabla u \| _ {H ^ {1} (\Omega)} \left\| h ^ {- \frac {1}{2}} [   [ v ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}}, \quad \forall v \in V _ {h} + H ^ {1} (\Omega).
$$

证明 我们有

$$
\begin{array}{l} \left| R _ {h} (u, v) \right| = \left| \langle [ v ] ], \{\kappa \nabla u - Q _ {h} (\kappa \nabla u) \} + \beta [ [ \kappa \nabla u - Q _ {h} (\kappa \nabla u) ] ] \rangle_ {\mathcal {E} _ {h} ^ {I D}} \right| \\ \lesssim \left(\sum_ {K \in \mathcal {M} _ {h}} h _ {K} \| \kappa \nabla u - Q _ {h} (\kappa \nabla u) \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {e \in \mathcal {E} _ {h} ^ {I D}} h _ {e} ^ {- 1} \| [ v ] \| _ {L ^ {2} (e)} ^ {2}\right) ^ {\frac {1}{2}} \\ \end{array}
$$

令 $\varphi_h|_K = \frac{1}{|K|}\int_K\kappa \nabla u\mathrm{d}\pmb{x}$ ，由局部迹不等式及有限元逆估计得

$$
\begin{array}{l} h _ {K} ^ {\frac {1}{2}} \left\| \kappa \nabla u - Q _ {h} (\kappa \nabla u) \right\| _ {\partial K} \\ \lesssim \| \kappa \nabla u - Q _ {h} (\kappa \nabla u) \| _ {K} + h _ {K} | \kappa \nabla u - Q _ {h} (\kappa \nabla u) | _ {H ^ {1} (K)} \\ \lesssim \| \kappa \nabla u - Q _ {h} (\kappa \nabla u) \| _ {K} + h _ {K} | \kappa \nabla u - \varphi_ {h} | _ {H ^ {1} (K)} + h _ {K} | \varphi_ {h} - Q _ {h} (\kappa \nabla u) | _ {H ^ {1} (K)} \\ \lesssim \| \kappa \nabla u - Q _ {h} (\kappa \nabla u) \| _ {K} + h _ {K} | \kappa \nabla u - \varphi_ {h} | _ {H ^ {1} (K)} + \| \varphi_ {h} - Q _ {h} (\kappa \nabla u) \| _ {K} \\ \lesssim \| \kappa \nabla u - Q _ {h} (\kappa \nabla u) \| _ {K} + h _ {K} | \kappa \nabla u - \varphi_ {h} | _ {H ^ {1} (K)} + \| \kappa \nabla u - \varphi_ {h} \| _ {K} \\ \lesssim \| \kappa \nabla u - \varphi_ {h} \| _ {K} + h _ {K} | \kappa \nabla u - \varphi_ {h} | _ {H ^ {1} (K)} \\ \lesssim h _ {K} | \kappa \nabla u | _ {H ^ {1} (K)}. \\ \end{array}
$$

将上面两个估计合并即得证明

注8.8如果 $\kappa$ 是分片线性的，那么在上面证明中取 $\varphi_h = \kappa \nabla I_h u$ 知 $\| \nabla u\|_{H^1 (\Omega)}$ 可以换为 $|u|_{H^2 (\Omega)}$

Strang 引理

由于LDG法不满足Galerkin正交性，所以Céa引理不成立，但利用引理8.7我可以证明如下引理：

引理8.9 假设 $|\beta_{e}| \lesssim 1, \gamma \gtrsim 1$ . $u$ 和 $u_{h}$ 分别是椭圆型问题(8.1)及其LDG离散(8.20)—(8.24)的解. 则

$$
\| u - u _ {h} \| _ {h} \lesssim \inf  _ {v _ {h} \in V _ {h}} \| u - v _ {h} \| _ {h} + \sup  _ {0 \neq w _ {h} \in V _ {h}} \frac {\left| B _ {h} (u , w _ {h}) - (f , w _ {h}) \right|}{\| | w _ {h} \| | _ {h}}.
$$

证明 由引理8.7及(8.32)，对任意 $v_{h}\in V_{h}$ 有

$$
\begin{array}{l} \| u _ {h} - v _ {h} \| _ {h} ^ {2} \lesssim B _ {h} \left(u _ {h} - v _ {h}, u _ {h} - v _ {h}\right) \\ = B _ {h} \left(u - v _ {h}, u _ {h} - v _ {h}\right) + \left(f, u _ {h} - v _ {h}\right) - B _ {h} \left(u, u _ {h} - v _ {h}\right) \\ \lesssim \| u - v _ {h} \| _ {h} \| u _ {h} - v _ {h} \| _ {h} + \\ \| u_{h} - v_{h}\|_{h}\sup_{0\neq w_{h}\in V_{h}}\frac{\left|B_{h}\left(u,w_{h}\right) - (f,w_{h})\right|}{\|w_{h}\|_{h}}. \\ \end{array}
$$

两边消去 $\| u_h - v_h\| _h$ ，再利用三角不等式即得证明

$H^{1}$ 误差估计

下面定理给出LDG法的 $H^{1}$ 误差估计

定理8.5 假设 $|\beta_{e}| \lesssim 1, \gamma \gtrsim 1, u$ 和 $u_{h}$ 分别是椭圆型问题(8.1）及其LDG离散(8.20)—(8.24）的解，则

$$
\| u - u _ {h} \| _ {h} \lesssim (1 + \gamma) ^ {\frac {1}{2}} h \| \nabla u \| _ {H ^ {1} (\Omega)}.
$$

证明 类似于引理8.3可得

$$
\| \| u - I _ {h} u \| \| _ {h} \lesssim (1 + \gamma) ^ {\frac {1}{2}} h | u | _ {H ^ {2} (\Omega)}.
$$

另外，由引理8.8知

$$
\left| B _ {h} \left(u - u _ {h}, w _ {h}\right) \right| \lesssim h \| \nabla u \| _ {H ^ {1} (\Omega)} \| w _ {h} \| _ {h}.
$$

将上面两个估计代入引理8.9即得证明

注8.9 对于协调的三角剖分，误差估计中的因子 $(1 + \gamma)^{\frac{1}{2}}$ 可以去掉

$L^2$ 误差估计

然后, 我们考虑 $L^2$ 误差估计. 设 $\boldsymbol{w}$ 为对偶问题 (8.14) 的解, 则

$$
\begin{array}{l} \left\| u - u _ {h} \right\| _ {L ^ {2} (\Omega)} ^ {2} = B _ {h} (w, u - u _ {h}) - R _ {h} (w, u - u _ {h}) \\ = B _ {h} \left(w - I _ {h} w, u - u _ {h}\right) + B _ {h} \left(u - u _ {h}, I _ {h} w\right) - R _ {h} \left(w, u - u _ {h}\right) \\ = B _ {h} \left(w - I _ {h} w, u - u _ {h}\right) + R _ {h} \left(u, I _ {h} w\right) - R _ {h} \left(w, u - u _ {h}\right). \\ \end{array}
$$

由引理8.7—8.8得

$$
\begin{array}{l} B _ {h} \left(w - I _ {h} w, u - u _ {h}\right) \lesssim \| | u - u _ {h} | | | _ {h} \| | | w - I _ {h} w | | | _ {h} \\ \lesssim (1 + \gamma) ^ {\frac {1}{2}} h \| u - u _ {h} \| _ {h} \| w \| _ {H ^ {2} (\Omega)} \\ \lesssim (1 + \gamma) ^ {\frac {1}{2}} h \| u - u _ {h} \| _ {h} \| u - u _ {h} \| _ {L ^ {2} (\Omega)}, \\ R _ {h} (u, I _ {h} w) \lesssim h \| \nabla u \| _ {H ^ {1} (\Omega)} \left\| h ^ {- \frac {1}{2}} [   [ I _ {h} w - w ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} \\ \lesssim \gamma^ {- \frac {1}{2}} h \| \nabla u \| _ {H ^ {1} (\Omega)} \| | w - I _ {h} w \| | _ {h} \\ \lesssim h ^ {2} \| \nabla u \| _ {H ^ {1} (\Omega)} \| u - u _ {h} \| _ {L ^ {2} (\Omega)}, \\ - R _ {h} (w, u - u _ {h}) \lesssim h \| w \| _ {H ^ {2} (\Omega)} \left\| h ^ {- \frac {1}{2}} [   [ u - u _ {h} ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} \\ \lesssim \gamma^ {- \frac {1}{2}} h \| u - u _ {h} \| _ {h} \| u - u _ {h} \| _ {L ^ {2} (\Omega)}. \\ \end{array}
$$

故

$$
\begin{array}{l} \| u - u _ {h} \| _ {L ^ {2} (\Omega)} ^ {2} \lesssim (1 + \gamma) ^ {\frac {1}{2}} h \| | u - u _ {h} | | | _ {h} \| u - u _ {h} \| _ {L ^ {2} (\Omega)} + \\ h ^ {2} \| \nabla u \| _ {H ^ {1} (\Omega)} \| u - u _ {h} \| _ {L ^ {2} (\Omega)}, \\ \end{array}
$$

结合定理8.5可得如下 $L^2$ 误差估计：

定理8.6 假设 $|\beta_{e}| \lesssim 1, \gamma \gtrsim 1, u$ 和 $u_{h}$ 分别是椭圆型问题(8.1）及其LDG离散(8.20)—(8.24）的解，则

$$
\| u - u _ {h} \| _ {L ^ {2} (\Omega)} \lesssim (1 + \gamma) h ^ {2} \| \nabla u \| _ {H ^ {1} (\Omega)}.
$$

注8.10 (a)LDG法与IPDG法相比， $u_{h}$ 的收敛阶相同，存储量前者较大LDG关于 $u_{h}$ 是对称格式且不要求加罚参数足够大；SIPG是对称格式，但要求加罚参数足够大；NIPG不要求加罚参数足够大，但格式不对称.

(b) $\sigma_h$ 与 $\kappa \nabla u_h$ 都是对 $\kappa \nabla u$ 的逼近. 由 (8.29) 及引理 8.6,

$$
\begin{array}{l} \left| \left\| \delta - \delta_ {h} \right\| - \left\| \nabla u - \nabla_ {h} u _ {h} \right\| \right| \leqslant \left\| L _ {h} u _ {h} \right\| \leqslant C _ {L} \left\| h ^ {- \frac {1}{2}} [   [ u _ {h} ]   ] \right\| _ {\mathcal {E} _ {h} ^ {I D}} \\ \leqslant C _ {L} \gamma^ {- \frac {1}{2}} \| u - u _ {h} \| _ {h}. \\ \end{array}
$$

可以看出，当 $\gamma$ 足够大的时候， $\sigma_h$ 与 $\kappa \nabla u_h$ 对 $\kappa \nabla u$ 的逼近是同阶的。由于 $\sigma_h$ 与 $u_h$ 取同阶的分片多项式，所以 $\sigma_h$ 的收敛阶一般不是最优的，但在某些特殊网格上， $\sigma_h$ 的逼近精度更高，甚至当 $\gamma$ 小的时候，可能达到最优阶。

# 8.3 杂交间断 Galerkin (HDG) 法

本节介绍HDG法[37].LDG法和IPDG法同连续有限元方法相比，有限元空间的定义更加灵活，但总自由度数一般要更多.HDG法的自由度定义在每个单元和单元的边上，但是每个单元上的自由度可以很容易地用其边上的自由度表示，所以最后的离散方程组仅仅是关于边上自由度的，有效地减少了总自由度数.所以，HDG法既具有DG法的灵活性，又具有总自由度数少的优点

# 8.3.1 离散格式

与前一节类似，引入中间变量 $\sigma = \kappa \nabla u$ 将椭圆型问题(8.1)改写为如下一阶偏微分方程组：

$$
\kappa^ {- 1} \sigma = \nabla u, - \nabla \cdot \sigma = f, \quad \boldsymbol {x} \in \Omega , \quad u | _ {\Gamma_ {D}} = 0, \quad (\sigma \cdot n) | _ {\Gamma_ {N}} = 0. \tag {8.37}
$$

用检验函数 $\tau$ 和 $v$ 分别乘前两个方程，并在任一单元 $K\in \mathcal{M}_h$ 上积分得

$$
\int_ {K} \kappa^ {- 1} \sigma \cdot \tau \mathrm {d} \boldsymbol {x} = - \int_ {K} u \nabla \cdot \tau \mathrm {d} \boldsymbol {x} + \int_ {\partial K} u \tau \cdot n _ {K} \mathrm {d} s, \tag {8.38}
$$

$$
\int_ {K} \sigma \cdot \nabla v \mathrm {d} x = \int_ {K} f v \mathrm {d} x + \int_ {\partial K} \sigma \cdot n _ {K} v \mathrm {d} s. \tag {8.39}
$$

为了离散化上面的变分形式，引入如下的逼近空间：

$$
V _ {h} := \left\{v _ {h} \in L ^ {2} (\Omega): v _ {h} | _ {K} \in P _ {1} (K), \forall K \in \mathcal {M} _ {h} \right\},
$$

$$
\Sigma_ {h} := \left(V _ {h}\right) ^ {d},
$$

$$
S _ {h} := \left\{v _ {h}: v _ {h} | _ {e} \in P _ {1} (e), \forall e \in \mathcal {E} _ {h} ^ {I N}, v _ {h} | _ {e} = 0, \forall e \in \mathcal {E} _ {h} ^ {D} \right\}.
$$

则求解椭圆型问题 (8.1) 的 HDG 法为: 求 $u_h \in V_h$ , $\sigma_h \in \Sigma_h$ 或求 $\hat{u}_h \in S_h$ 使得下面 (8.40) ——(8.43) 成立

$$
\int_ {K} \kappa^ {- 1} \sigma_ {h} \cdot \tau_ {h} d x = - \int_ {K} u _ {h} \nabla \cdot \tau_ {h} d x + \int_ {\partial K} \hat {u} _ {h} \tau_ {h} \cdot n _ {K} d s, \quad \forall \tau_ {h} \in \Sigma_ {h}, \tag {8.40}
$$

$$
\int_ {K} \sigma_ {h} \cdot \nabla v _ {h} d x = \int_ {K} f v _ {h} d x + \int_ {\partial K} \hat {\sigma} _ {h} \cdot n _ {K} v _ {h} d s, \quad \forall v _ {h} \in V _ {h}, K \in \mathcal {M} _ {h}. \tag {8.41}
$$

其中 $\hat{\sigma}_h$ 和 $\hat{u}_h$ 称为数值迹，分别是 $\kappa \nabla u$ 和 $u$ 的近似，如下定义：

$$
\hat {\sigma} _ {h} \mid_ {\partial K} := \sigma_ {h} + \alpha (\hat {u} _ {h} - u _ {h}) n _ {K}, \quad \forall K \in \mathcal {M} _ {h}. \tag {8.42}
$$

并且要求

$$
\left. \llbracket \hat {\sigma} _ {h} \rrbracket \right| _ {e} = 0, \quad \forall e \in \mathcal {E} _ {h} ^ {I N}, \tag {8.43}
$$

其中 $\alpha > 0$ 为稳定化参数或称为加罚参数

注意到 (8.40) (8.41) 与 LDG 格式中的 (8.20) — (8.22) 基本相同. 而且 (8.42) (8.43) 也可以改成类似 (8.23) (8.24) 的形式. 事实上, 记

$$
\alpha_ {*} = \left\{ \begin{array}{l l} \alpha , & e \in \mathcal {E} _ {h} ^ {I}, \\ 2 \alpha , & e \in \mathcal {E} _ {h} ^ {D}, \\ \frac {1}{2} \alpha , & e \in \mathcal {E} _ {h} ^ {N}. \end{array} \right. \tag {8.44}
$$

由(8.42)得

$$
\llbracket \hat {\sigma} _ {h} \rrbracket = \llbracket \sigma_ {h} \rrbracket + 2 \alpha_ {*} (\hat {u} _ {h} - \{u _ {h} \}), \quad e \in \mathcal {E} _ {h} ^ {I N}.
$$

再由 (8.43) 得

$$
\hat {u} _ {h} = \left\{u _ {h} \right\} - \frac {1}{2 \alpha_ {*}} \llbracket \sigma_ {h} \rrbracket , \quad \forall e \in \mathcal {E} _ {h} ^ {I N}. \tag {8.45}
$$

另外，再由 (8.42) 可得

$$
\hat {\sigma} _ {h} = \left\{\sigma_ {h} \right\} - \frac {\alpha_ {*}}{2} \llbracket u _ {h} \rrbracket , \quad \forall e \in \mathcal {E} _ {h} ^ {I D}. \tag {8.46}
$$

# 8.3.2 变分形式I

为了理论分析的方便, 本小节将HDG格式(8.40)一(8.43)改写为关于 $u_{h}, \sigma_{h}$ 的变分公式. 首先, 将(8.40)关于 $K \in \mathcal{M}_{h}$ 求和并利用(8.25)得

$$
\left(\kappa^ {- 1} \sigma_ {h}, \tau_ {h}\right) = - \left(u _ {h}, \nabla_ {h} \cdot \tau_ {h}\right) + \langle \hat {u} _ {h}, [   [ \tau_ {h} ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I N}}. \tag {8.47}
$$

将(8.45)代入(8.47)整理得

$$
\left(\kappa^ {- 1} \sigma_ {h}, \tau_ {h}\right) + \left(u _ {h}, \nabla_ {h} \cdot \tau_ {h}\right) + \left\langle \frac {1}{2 \alpha_ {*}} \llbracket \sigma_ {h} \rrbracket - \left\{u _ {h} \right\}, \llbracket \tau_ {h} \rrbracket \right\rangle_ {\mathcal {E} _ {h} ^ {I N}} = 0. \tag {8.48}
$$

同理, 对 (8.41) 关于 $K \in \mathcal{M}_h$ 求和并利用 (8.25) 及得

$$
\left(\sigma_ {h}, \nabla_ {h} v _ {h}\right) = (f, v _ {h}) + \langle \{\hat {\sigma} _ {h} \}, [   [ v _ {h} ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I D}}.
$$

由(8.46)得

$$
\left. \left(\sigma_ {h}, \nabla_ {h} v _ {h}\right) + \left\langle \frac {\alpha_ {*}}{2} \llbracket u _ {h} \rrbracket - \left\{\sigma_ {h} \right\}, \llbracket v _ {h} \rrbracket \right\rangle_ {\mathcal {E} _ {h} ^ {I D}} = (f, v _ {h}) \right.. \tag {8.49}
$$

对 $u,v\in V_h + H^1 (\varOmega),\sigma ,\tau \in \Sigma_h + \left(H^1 (\varOmega)\right)^d$ ，定义

$$
\begin{array}{l} A (u, \sigma ; v, \tau) = (\kappa^ {- 1} \sigma , \tau) + (u, \nabla_ {h} \cdot \tau) + (\sigma , \nabla_ {h} v) + \\ \left\langle \frac {\alpha_ {*}}{2} [   [ u ]   ] - \{\sigma \}, [   [ v ]   ] \right\rangle_ {\mathcal {E} _ {h} ^ {I D}} + \left\langle \frac {1}{2 \alpha_ {*}} [   [ \sigma ]   ] - \{u \}, [   [ \tau ]   ] \right\rangle_ {\mathcal {E} _ {h} ^ {I N}}. \tag {8.50} \\ \end{array}
$$

将(8.48)和(8.49)相加知，得到HDG法关于解 $u_{h},\sigma_{h}$ 的变分形式：

$$
A \left(u _ {h}, \sigma_ {h}; v _ {h}, \tau_ {h}\right) = (f, v _ {h}), \quad \forall v _ {h} \in V _ {h}, \tau_ {h} \in \Sigma_ {h}. \tag {8.51}
$$

为了理论分析的方便，下面给出 $A(u,\sigma ;v,\tau)$ 的另外两种表示.由（8.25）易知

$$
(u, \nabla_ {h} \cdot \tau) = - (\nabla_ {h} u, \tau) + \langle [   [ u ]   ], \{\tau \} \rangle_ {\mathcal {E} _ {h} ^ {I D}} + \langle \{u \}, [   [ \tau ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I N}},
$$

$$
(\sigma , \nabla_ {h} v) = - (\nabla_ {h} \cdot \sigma , v) + \langle \{\sigma \}, [   [ v ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I D}} + \langle [   [ \sigma ]   ], \{v \} \rangle_ {\mathcal {E} _ {h} ^ {I N}}.
$$

代入到 (8.50) 得

$$
\begin{array}{l} A (u, \sigma ; v, \tau) = (\kappa^ {- 1} \sigma , \tau) - (\nabla_ {h} u, \tau) + (\sigma , \nabla_ {h} v) - \\ \langle \{\sigma \}, [   [ v ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I D}} + \langle [   [ u ]   ], \{\tau \} \rangle_ {\mathcal {E} _ {h} ^ {I D}} + \\ \left\langle \frac {\alpha_ {*}}{2} [   [ u ]   ], [   [ v ]   ] \right\rangle_ {\mathcal {E} _ {h} ^ {I D}} + \left\langle \frac {1}{2 \alpha_ {*}} [   [ \sigma ]   ], [   [ \tau ]   ] \right\rangle_ {\mathcal {E} _ {h} ^ {I N}} \tag {8.52} \\ \end{array}
$$

及

$$
\begin{array}{l} A (u, \sigma ; v, \tau) = (\kappa^ {- 1} \sigma , \tau) - (\nabla_ {h} u, \tau) - (\nabla_ {h} \cdot \sigma , v) + \\ \left\langle \llbracket u \rrbracket , \frac {\alpha_ {*}}{2} \llbracket v \rrbracket + \{\tau \} \right\rangle_ {\mathcal {E} _ {h} ^ {I D}} + \left\langle \llbracket \sigma \rrbracket , \frac {1}{2 \alpha_ {*}} \llbracket \tau \rrbracket + \{v \} \right\rangle_ {\mathcal {E} _ {h} ^ {I N}}. \tag {8.53} \\ \end{array}
$$

定理8.7 设 $\alpha > 0$ ，则HDG法(8.40)一(8.43)存在唯一解

证明 只需证明 $f = 0$ 时必有 $u_{h} \equiv 0, \sigma_{h} \equiv 0, \hat{u}_{h} \equiv 0$ . 在 (8.51) 中取 $f = 0, v_{h} = u_{h}, \tau_{h} = \sigma_{h}$ , 并利用 (8.52) 得

$$
\left\| \kappa^ {- \frac {1}{2}} \sigma_ {h} \right\| _ {L ^ {2} (\Omega)} ^ {2} + \left\langle \frac {\alpha_ {*}}{2} [   [ u _ {h} ]   ], [   [ u _ {h} ]   ] \right\rangle_ {\mathcal {E} _ {h} ^ {I D}} + \left\langle \frac {1}{2 \alpha_ {*}} [   [ \sigma_ {h} ]   ], [   [ \sigma_ {h} ]   ] \right\rangle_ {\mathcal {E} _ {h} ^ {I N}} = 0.
$$

因此

$$
\sigma_ {h} \equiv 0; \quad u _ {h} | _ {\Gamma_ {D}} = 0; \quad u _ {h} \in C (\bar {\Omega}).
$$

代入到 (8.51) 并利用 (8.53) 知, $(\nabla_h u_h, \tau_h) = 0$ , $\forall \tau_h \in \Sigma_h$ , 从而 $\nabla_h u_h = 0$ , 故 $u_h$ 在 $\mathcal{M}_h$ 上是分片常数, 只能有 $u_h \equiv 0$ , 进而知 $\hat{u}_h \equiv 0$ (参见 (8.45)). 证毕. □

# 8.3.3 误差估计

首先，由(8.37)及(8.53)易知，椭圆型问题的解 $u,\sigma$ 满足：

$$
A (u, \sigma ; v _ {h}, \tau_ {h}) = (f, v _ {h}), \quad \forall v _ {h} \in V _ {h}, \tau_ {h} \in \Sigma_ {h}. \tag {8.54}
$$

即HDG法的变分形式(8.51)与椭圆型问题(8.37）是相容的，从而得Galerkin正交性：

$$
A \left(u - u _ {h}, \sigma - \sigma_ {h}; v _ {h}, \tau_ {h}\right) = 0, \quad \forall v _ {h} \in V _ {h}, \tau_ {h} \in \Sigma_ {h}. \tag {8.55}
$$

对于给定常数 $\beta$ ，我们引入到 $V_{h}\times \Sigma_{h}$ 的投影算子 $\Pi_h^\beta$ .定义 $\Pi_h^\beta (u,\sigma) = \left(\Pi_1^\beta u,\Pi_2^\beta \sigma\right)\in$ $V_{h}\times \Sigma_{h}$ 满足： $\forall K\in \mathcal{M}_h$ ， $e\subset \partial K$ ， $e\in \mathcal{E}_h$

$$
\begin{array}{l} \left(\Pi_ {1} ^ {\beta} u, v\right) _ {K} = (u, v) _ {K}, \quad \forall v \in P _ {0} (K), (8.56) \\ \left(\Pi_ {2} ^ {\beta} \sigma , \tau\right) _ {K} = (\sigma , \tau) _ {K}, \quad \forall \tau \in P _ {0} (K) ^ {d}, (8.57) \\ \left\langle \Pi_ {2} ^ {\beta} \sigma \cdot n _ {K} - \beta \Pi_ {1} ^ {\beta} u, \mu \right\rangle_ {e} = \left\langle \sigma \cdot n _ {K} - \beta u, \mu \right\rangle_ {e}, \quad \forall \mu \in P _ {1} (e). (8.58) \\ \end{array}
$$

引理8.10 设 $1 \leqslant s, t \leqslant 2$ ，则 $\forall K \in \mathcal{M}_h$

$$
\begin{array}{l} | \beta | \| \Pi_ {1} ^ {\beta} u - u \| _ {K} \lesssim | \beta | h _ {K} ^ {s} | u | _ {H ^ {s} (K)} + h _ {K} ^ {t} | \nabla \cdot \sigma | _ {H ^ {t - 1} (K)}, (8.59) \\ \left\| \Pi_ {2} ^ {\beta} \sigma - \sigma \right\| _ {K} \lesssim | \beta | h _ {K} ^ {s} | u | _ {H ^ {s} (K)} + h _ {K} ^ {t} | \sigma | _ {H ^ {t} (K)}. (8.60) \\ \end{array}
$$

证明 分为以下几步证明：

(i) 先证 $\Pi_h^\beta$ 的存在唯一性. 设 $u = 0$ , $\sigma = 0$ . 在 (8.56) 和 (8.57) 中分别取 $v = \nabla \cdot \Pi_2^\beta \sigma$ , $\tau = \nabla \Pi_1^\beta u$ 并相加得:

$$
\left\langle \Pi_ {2} ^ {\beta} \sigma \cdot n _ {K}, \Pi_ {1} ^ {\beta} u \right\rangle_ {\partial K} = \left(\Pi_ {1} ^ {\beta} u, \nabla \cdot \Pi_ {2} ^ {\beta} \sigma\right) _ {K} + \left(\Pi_ {2} ^ {\beta} \sigma , \nabla \Pi_ {1} ^ {\beta} u\right) _ {K} = 0.
$$

由(8.58)得：在 $\partial K$ 上 $\Pi_2^\beta \sigma \cdot n_K = \beta \Pi_1^\beta u.$ 代入上式得 $\left(\Pi_1^\beta u\right)\big|_{\partial K} = 0,\left(\Pi_2^\beta \sigma \cdot n_K\right)\big|_{\partial K} = 0,$ 推出 $\Pi_1^\beta u$ 和 $\Pi_2^\beta \sigma$ 在单元 $K$ 的每个顶点处为零，从而 $\left(\Pi_1^\beta u\right)\big|_K = 0,\left(\Pi_2^\beta \sigma\right)_K = 0.$ 存在唯一性得证

(ii) 由 (i) 的证明可知: 如果 $u \in V_h$ , $\sigma \in \Sigma_h$ , 则 $\Pi_1^\beta u = u$ , $\Pi_2^\beta \sigma = \sigma$ .  
(iii) 稳定性. 首先, 由尺度变换技巧易证

$$
\| v _ {h} \| _ {K} \approx h _ {K} ^ {\frac {1}{2}} \| v _ {h} \| _ {\partial K}, \forall v _ {h} \in V _ {h}; \| \tau_ {h} \| _ {K} \approx h _ {K} ^ {\frac {1}{2}} \| \tau_ {h} \cdot n _ {K} \| _ {\partial K}, \forall \tau_ {h} \in \Sigma_ {h}.
$$

其中, “ $A \approx B$ ” 表示存在不依赖步长 $h$ 的正常数 $C_1$ 和 $C_2$ , 使得 $A \leqslant C_1B$ 且 $B \leqslant C_2A$ . 设 $Q_K: L^2(K) \mapsto P_1(K)$ 为单元 $K$ 上的 $L^2$ 正交投影算子. 在 (8.58) 中取 $\mu = Q_K u - \Pi_1^\beta u \in P_1(K)$ 并注意到 $\int_K \mu = 0$ 得:

$$
\begin{array}{l} \beta \| \Pi_ {1} ^ {\beta} u \| _ {\partial K} ^ {2} = \beta \left\langle \Pi_ {1} ^ {\beta} u, Q _ {K} u \right\rangle_ {\partial K} - \beta \langle u, \mu \rangle_ {\partial K} + \left\langle \left(\sigma - \Pi_ {2} ^ {\beta} \sigma\right) \cdot n _ {K}, \mu \right\rangle_ {\partial K} \\ = \beta \left\langle \Pi_ {1} ^ {\beta} u, Q _ {K} u \right\rangle_ {\partial K} - \beta \langle u, \mu \rangle_ {\partial K} + \left(\nabla \cdot \left(\sigma - \Pi_ {2} ^ {\beta} \sigma\right), \mu\right) _ {K} \\ = \beta \left\langle \Pi_ {1} ^ {\beta} u, Q _ {K} u \right\rangle_ {\partial K} - \beta \langle u, \mu \rangle_ {\partial K} + (\nabla \cdot \sigma , \mu) _ {K} \\ \lesssim \beta h _ {K} ^ {- 1} \| \Pi_ {1} ^ {\beta} u \| _ {K} \| Q _ {K} u \| _ {K} + \beta h _ {K} ^ {- \frac {1}{2}} \| u \| _ {\partial K} \| \mu \| _ {K} + \| \nabla \cdot \sigma \| _ {K} \| \mu \| _ {K}. \\ \end{array}
$$

故

$$
\begin{array}{l} | \beta | ^ {2} \| \Pi_ {1} ^ {\beta} u \| _ {K} ^ {2} \lesssim | \beta | ^ {2} h _ {K} \| \Pi_ {1} ^ {\beta} u \| _ {\partial K} ^ {2} \lesssim | \beta | ^ {2} \| \Pi_ {1} ^ {\beta} u \| _ {K} \| Q _ {K} u \| _ {K} + \\ \left(| \beta | ^ {2} h _ {K} ^ {\frac {1}{2}} \| u \| _ {\partial K} + h _ {K} | \beta | \| \nabla \cdot \sigma \| _ {K}\right) \left(\| Q _ {K} u \| _ {K} + \left\| \Pi_ {1} ^ {\beta} u \right\| _ {K}\right). \\ \end{array}
$$

易得

$$
\left| \beta \right| ^ {2} \left\| \Pi_ {1} ^ {\beta} u \right\| _ {K} ^ {2} \lesssim \left| \beta \right| ^ {2} \left\| Q _ {K} u \right\| _ {K} ^ {2} + \left| \beta \right| ^ {2} h _ {K} \| u \| _ {\partial K} ^ {2} + h _ {K} ^ {2} \| \nabla \cdot \sigma \| _ {K} ^ {2},
$$

从而，由局部迹不等式可得 $\Pi_1^\beta u$ 的估计

$$
\left| \beta \right| \left\| \Pi_ {1} ^ {\beta} u \right\| _ {K} \lesssim \left| \beta \right| \left(\| u \| _ {K} + h _ {K} \| \nabla u \| _ {K}\right) + h _ {K} \| \nabla \cdot \sigma \| _ {K}.
$$

另外, 在 (8.58) 中取 $\mu = \Pi_2^\beta \sigma \cdot n_K - \beta \Pi_1^\beta u$ 得:

$$
\left| \left| \Pi_ {2} ^ {\beta} \sigma \cdot n _ {K} - \beta \Pi_ {1} ^ {\beta} u \right| \right| _ {\partial K} \leqslant \| \sigma \| _ {\partial K} + | \beta | \| u \| _ {\partial K}.
$$

从而可得 $\Pi_2^\beta \sigma$ 的估计

$$
\begin{array}{l} \left\| \Pi_ {2} ^ {\beta} \sigma \right\| _ {K} \lesssim h _ {K} ^ {\frac {1}{2}} \left\| \Pi_ {2} ^ {\beta} \sigma \cdot n _ {K} \right\| _ {\partial K} \\ \lesssim | \beta | \left(\| u \| _ {K} + h _ {K} \| \nabla u \| _ {K}\right) + \| \sigma \| _ {K} + h _ {K} | \sigma | _ {H ^ {1} (K)}. \\ \end{array}
$$

(iv) 由 (ii) 得 $(u, \sigma) - \Pi_h^\beta(u, \sigma) = (u, \sigma) - (v_h, \tau_h) + \Pi_h^\beta(v_h - u, \tau_h - \sigma), \forall v_h \in V_h, \tau_h \in \Sigma_h$ . 再由 (iii) 得

$$
\begin{array}{l} | \beta | \| u - \Pi_ {1} ^ {\beta} u \| _ {K} = | \beta | \| u - v _ {h} + \Pi_ {1} ^ {\beta} (v _ {h} - u) \| _ {K} \\ \lesssim | \beta | \left(\| u - v _ {h} \| _ {K} + h _ {K} \| \nabla (u - v _ {h}) \| _ {K}\right) + h _ {K} \| \nabla \cdot (\sigma - \tau_ {h}) \| _ {K}. \\ \end{array}
$$

为了证明(8.59)，可以如下选取 $v_{h}$ 和 $\tau_h$ ：对 $s = 1$ ，取 $v_{h}|_{K}$ 为 $u$ 的积分平均 $u_{K}$ ，对 $s = 2$ 取 $v_{h} = I_{h}u$ ；对 $t = 1$ ，取 $\tau_h|_K = 0$ ，对 $t = 2$ ，取 $\tau_h\in \Sigma_h$ ，使得 $(\nabla \cdot \tau_h)_K = (\nabla \cdot \sigma)_K.$ 类似可证(8.60).证毕. □

简记 $\Pi_{i} = \Pi_{i}^{\alpha}, i = 1,2.$ 将误差分解为

$$
u - u _ {h} = (u - \Pi_ {1} u) - (u _ {h} - \Pi_ {1} u) =: q - s _ {h};
$$

$$
\sigma - \sigma_ {h} = (\sigma - \Pi_ {2} \sigma) - (\sigma_ {h} - \Pi_ {2} \sigma) =: \rho - \theta_ {h}.
$$

由 $\Pi_h^\beta$ 的定义知：

$$
(q, v) _ {K} = 0, \forall v \in P _ {0} (K); \quad (\rho , \tau) _ {K} = 0, \forall \tau \in P _ {0} (K) ^ {d}, \quad K \in \mathcal {M} _ {h}; \tag {8.61}
$$

$$
\langle [ \rho ], \mu \rangle_ {e} = 2 \alpha_ {*} \langle \{q \}, \mu \rangle_ {e}, \quad \forall \mu \in P _ {1} (e), \quad e \in \mathcal {E} _ {h} ^ {I N}; \tag {8.62}
$$

$$
2 \langle \{\rho \}, \mu n \rangle_ {e} = \alpha_ {*} \langle [   [ q ]   ], \mu n \rangle_ {e}, \quad \forall \mu \in P _ {1} (e), \quad e \in \mathcal {E} _ {h} ^ {I D}. \tag {8.63}
$$

由(8.55)(8.50)和(8.61)—(8.63)得

$$
\begin{array}{l} A \left(s _ {h}, \theta_ {h}; v _ {h}, \tau_ {h}\right) = A (q, \rho ; v _ {h}, \tau_ {h}) = (\kappa^ {- 1} \rho , \tau_ {h}) + (q, \nabla_ {h} \cdot \tau_ {h}) + (\rho , \nabla_ {h} v _ {h}) \tag {8.64} \\ = \left(\kappa^ {- 1} \rho , \tau_ {h}\right), \quad \forall v _ {h} \in V _ {h}, \tau_ {h} \in \Sigma_ {h}. \\ \end{array}
$$

定义离散范数

$$
\| \left| v, \tau \right| \| _ {h} ^ {2} := A (v, \tau ; v, \tau) = (\kappa^ {- 1} \tau , \tau) + \left\langle \frac {\alpha_ {*}}{2} [   [ v ]   ], [   [ v ]   ] \right\rangle_ {\mathcal {E} _ {h} ^ {I D}} + \left\langle \frac {1}{2 \alpha_ {*}} [   [ \tau ]   ], [   [ \tau ]   ] \right\rangle_ {\mathcal {E} _ {h} ^ {I N}}. \tag {8.65}
$$

从而

$$
\| s _ {h}, \theta_ {h} \| _ {h} ^ {2} = A (s _ {h}, \theta_ {h}; s _ {h}, \theta_ {h}) = (\kappa^ {- 1} \rho , \theta_ {h}) \leqslant \| \kappa^ {- \frac {1}{2}} \rho \| \| \kappa^ {- \frac {1}{2}} \theta_ {h} \|.
$$

得

$$
\| s _ {h}, \theta_ {h} \| _ {h} \leqslant \| \kappa^ {- \frac {1}{2}} \rho \|. \tag {8.66}
$$

再利用引理8.10可得如下定理：

定理8.8 假设 $0 < \alpha \lesssim 1$ ，则

$$
\| \sigma - \sigma_ {h} \| \lesssim \| \rho \| \lesssim h ^ {2} \left(| u | _ {H ^ {2} (\Omega)} + | \sigma | _ {H ^ {2} (\Omega)}\right).
$$

注8.11 若 $0 < \alpha \lesssim h^{-1}$ , 则 $\| \sigma - \sigma_h \| \lesssim h |u|_{H^2(\Omega)}$ .

下面用对偶论证技巧推导 $u_{h}$ 的 $L^2$ 误差估计. 引入对偶问题:

$$
\varphi = - \kappa \nabla w, \quad \nabla \cdot \varphi = s _ {h}, \quad \boldsymbol {x} \in \Omega , \quad w | _ {\Gamma_ {D}} = 0, \quad (\varphi \cdot n) | _ {\Gamma_ {N}} = 0. \tag {8.67}
$$

需要注意的是，上面引入中间变量时的 $\nabla$ 算子前面的符号和原问题一阶组形式 $\nabla$ 前面的符号正好相反.我们有正则性估计：

$$
\| w \| _ {H ^ {2} (\Omega)} + \| \varphi \| _ {H ^ {1} (\Omega)} \lesssim \| s _ {h} \|. \tag {8.68}
$$

由 (8.50) 得

$$
A (v, \tau ; w, \varphi) = (v, \nabla \cdot \varphi) = (v, s _ {h}). \tag {8.69}
$$

简记 $\Pi_i^* = \Pi_i^{-\alpha}, i = 1,2, q^* = w - \Pi_1^* w, \rho^* = \varphi - \Pi_2^* \varphi.$ 类似 (8.62) (8.63) 有

$$
\begin{array}{l} \langle [   [ \rho^ {*} ]   ], \mu \rangle_ {e} = - 2 \alpha_ {*} \langle \{q ^ {*} \}, \mu \rangle_ {e}, \quad \forall \mu \in P _ {1} (e), e \in \mathcal {E} _ {h} ^ {I N}; \\ 2 \left\langle \left\{\rho^ {*} \right\}, \mu n \right\rangle_ {e} = - \alpha_ {*} \left\langle \llbracket q ^ {*} \rrbracket , \mu n \right\rangle_ {e}, \quad \forall \mu \in P _ {1} (e), e \in \mathcal {E} _ {h} ^ {I D}. \\ \end{array}
$$

在 (8.69) 中取 $v = s_h, \tau = \theta_h$ 并利用 (8.64) (8.53), 上面两式及 (8.56) (8.57) 得:

$$
\begin{array}{l} \left\| s _ {h} \right\| ^ {2} = A \left(s _ {h}, \theta_ {h}; w, \varphi\right) = A \left(s _ {h}, \theta_ {h}; \Pi_ {1} ^ {*} w, \Pi_ {2} ^ {*} \varphi\right) + A \left(s _ {h}, \theta_ {h}; q ^ {*}, \rho^ {*}\right) \\ = \left(\kappa^ {- 1} \rho , \Pi_ {2} ^ {*} \varphi\right) + \left(\kappa^ {- 1} \theta_ {h}, \rho^ {*}\right) \\ = \left(\kappa^ {- 1} (\rho - \theta_ {h}), \Pi_ {2} ^ {*} \varphi - \varphi\right) - \left(\rho , \nabla_ {h} (w - I _ {h} w)\right). \\ \end{array}
$$

由 (8.66) (8.68), 引理 8.10, 当 $\alpha h \lesssim 1$ 时,

$$
\begin{array}{l} \left\| s _ {h} \right\| ^ {2} \lesssim \| \rho \| \left(\| \varphi - \Pi_ {2} ^ {*} \varphi \| + h | w | _ {H ^ {2} (\Omega)}\right) \\ \lesssim \| \rho \| \left(\left(| \alpha | h ^ {2} + h\right) | w | _ {H ^ {2} (\Omega)} + h | \varphi | _ {H ^ {1} (\Omega)}\right) \lesssim h \| s _ {h} \| \| \rho \|. \\ \end{array}
$$

得如下定理：

定理8.9 假设 $0 < \alpha \lesssim 1$ ，则

$$
\| \Pi_ {1} u - u _ {h} \| \lesssim h \| \rho \| \lesssim h ^ {3} \left(| u | _ {H ^ {2} (\Omega)} + | \sigma | _ {H ^ {2} (\Omega)}\right).
$$

进一步，若 $\alpha \sim 1$ ，则

$$
\left\| u - u _ {h} \right\| \lesssim \| q \| + h \| \rho \| \lesssim h ^ {2} \left(\left| u \right| _ {H ^ {2} (\Omega)} + \left| \sigma \right| _ {H ^ {1} (\Omega)} + \left| \nabla \cdot \sigma \right| _ {H ^ {1} (\Omega)}\right).
$$

注8.12 (a) 注意到 $u_h$ 与 $\Pi_1 u$ 之间的误差有超收敛。我们可以利用这一点，通过后处理来提高 $u_h$ 的精度。

(b) 若 $\alpha \approx h^{-1}$ , 则 $\|u - u_h\| \lesssim h^2(|u|_{H^2(\Omega)} + |\sigma|_{H^1(\Omega)})$

# 8.3.4 后处理

本小节利用较少的计算量对 $u_{h}$ 进行后处理，得到精确解 $\pmb{u}$ 的分片二次近似，提高收敛阶.再由 $\int_{\Omega}|\nabla u|^{2}d\Omega = 0,\int_{\Omega}|\nabla u|^{2}d\Omega = 0,\int_{\Omega}|\nabla u|^{2}d\Omega = 0$

记 $m_K(v) = \frac{1}{|K|}\int_K v$ 为 $v$ 在单元 $K$ 上的积分平均. 显然 $m_K\pi_1u = m_Ku,$ 从而 $u_h$ 在 $\mathcal{M}_h$ 上的分片积分平均函数对精确解 $u$ 的分片积分平均函数具有超逼近性质. 记 $P_2^0 (K)\coloneqq \{p\in P_2(K):m_K(p) = 0\} .$ 定义 $u_{h}^{*}$ 如下： $u_{h}^{*}|_{K}\in P_{2}(K),$

$$
\left\{ \begin{array}{l l} \left(\kappa \nabla u _ {h} ^ {*}, \nabla v\right) _ {K} = (f, v) _ {K} + \langle \hat {\sigma} _ {h} \cdot n _ {K}, v \rangle_ {\partial K}, & \forall v \in P _ {2} ^ {0} (K), \\ m _ {K} \left(u _ {h} ^ {*}\right) = m _ {K} \left(u _ {h}\right), & \forall K \in \mathcal {M} _ {h}. \end{array} \right. \tag {8.70}
$$

记 $V_{h,2} \coloneqq \prod_{K \in \mathcal{M}_h} P_2(K)$ . 下面定理给出 $u_h^*$ 的高阶估计.

定理8.10 假设 $\alpha \sim 1$ ，则

$$
\left\| u - u _ {h} ^ {*} \right\| \lesssim h ^ {3} \left(\left| u \right| _ {H ^ {2} (\Omega)} + \left| \sigma \right| _ {H ^ {2} (\Omega)}\right).
$$

证明 $\forall v_{h}\in V_{h,2}$ ，取 $v = (I - m_K)(v_h - u_h^*)$ 得：

$$
\begin{array}{l} \left\| \kappa^ {\frac {1}{2}} \nabla v \right\| _ {K} ^ {2} = \left(\kappa \nabla \left(v _ {h} - u _ {h} ^ {*}\right), \nabla v\right) _ {K} \\ = (\kappa \nabla (v _ {h} - u), \nabla v) _ {K} + \langle (\sigma - \sigma_ {h}) \cdot n _ {K}, v \rangle_ {\partial K} + \langle (\sigma_ {h} - \hat {\sigma} _ {h}) \cdot n _ {K}, v \rangle_ {\partial K} \\ = (\kappa \nabla (v _ {h} - u), \nabla v) _ {K} + (\nabla \cdot (\sigma - \sigma_ {h}), v) _ {K} + (\sigma - \sigma_ {h}, \nabla v) _ {K} + \\ \left\langle \left(\sigma_ {h} - \hat {\sigma} _ {h}\right) \cdot n _ {K}, v \right\rangle_ {\partial K} \\ = (\kappa \nabla (v _ {h} - u), \nabla v) _ {K} + (\nabla \cdot (\sigma - \tau_ {h}), v) _ {K} + (\sigma - \sigma_ {h}, \nabla v) _ {K} + \\ \left\langle \left(\sigma_ {h} - \hat {\sigma} _ {h}\right) \cdot n _ {K}, v \right\rangle_ {\partial K}, \\ \end{array}
$$

其中 $\tau_h\in \Sigma_h$ ，这里我们用到了 $\int_{K}v = 0.$ 从而由 $\| v\| _K\lesssim h_K\| \nabla v\| _K$ 及 $\| v\|_{\partial K}\lesssim$ $h_K^{\frac{1}{2}}\| \nabla v\| _K$ 得

$$
\begin{array}{l} h _ {K} ^ {- 1} \| v \| _ {K} \lesssim \| \nabla v \| _ {K} \lesssim \| \nabla (u - v _ {h}) \| _ {K} + h _ {K} \| \nabla \cdot (\sigma - \tau_ {h}) \| _ {K} + \| \sigma - \sigma_ {h} \| _ {K} + \\ h _ {K} ^ {\frac {1}{2}} \left\| \left(\sigma_ {h} - \hat {\sigma} _ {h}\right) \cdot n _ {K} \right\| _ {\partial K}. \\ \end{array}
$$

从而由三角不等式得

$$
\begin{array}{l} \left\| u - u _ {h} ^ {*} \right\| _ {K} \leqslant \left\| u - v _ {h} \right\| _ {K} + \left\| v \right\| _ {K} + \left\| m _ {K} \left(v _ {h} - u _ {h} ^ {*}\right) \right\| _ {K} \\ \lesssim \| u - v _ {h} \| _ {K} + h _ {K} \| \nabla (u - v _ {h}) \| _ {K} + h _ {K} ^ {2} \| \nabla \cdot (\sigma - \tau_ {h}) \| _ {K} + \\ \left\| m _ {K} \left(v _ {h} - u _ {h} ^ {*}\right) \right\| _ {K} + h _ {K} \left\| \sigma - \sigma_ {h} \right\| _ {K} + h _ {K} ^ {\frac {3}{2}} \left\| \left(\sigma_ {h} - \hat {\sigma} _ {h}\right) \cdot n _ {K} \right\| _ {\partial K}. \\ \end{array}
$$

最后，注意到

$$
\left\| m _ {K} \left(v _ {h} - u _ {h} ^ {*}\right) \right\| _ {K} = \left\| m _ {K} \left(v _ {h} - u\right) + m _ {K} \left(\Pi_ {1} u - u _ {h}\right) \right\| _ {K} \leqslant \left\| u - v _ {h} \right\| _ {K} + \left\| \Pi_ {1} u - u _ {h} \right\| _ {K}
$$

及（由（8.46）（8.66）及局部迹不等式）

$$
\begin{array}{l} \left(\sum_ {K \in \mathcal {M} _ {h}} \| (\sigma_ {h} - \hat {\sigma} _ {h}) \cdot n _ {K} \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \lesssim \| [ \sigma_ {h} ] \| _ {\mathcal {E} _ {h} ^ {I N}} + \alpha \| [ [ u _ {h} ] ] \| _ {\mathcal {E} _ {h} ^ {I D}} \\ = \| \llbracket \sigma_ {h} - \tau_ {h} + \tau_ {h} - \sigma \rrbracket \| _ {\mathcal {E} _ {h} ^ {I N}} + \alpha \| \llbracket u _ {h} - v _ {h} + v _ {h} - u \rrbracket \| _ {\mathcal {E} _ {h} ^ {I D}} \\ \lesssim h ^ {- \frac {1}{2}} \left(\| \sigma_ {h} - \tau_ {h} \| + \| \sigma - \tau_ {h} \| + \| u _ {h} - v _ {h} \| + \| u - v _ {h} \|\right) + \\ h ^ {\frac {1}{2}} \left(\left\| \nabla_ {h} (\sigma - \tau_ {h}) \right\| + \left\| \nabla_ {h} (u - v _ {h}) \right\|\right) \\ \lesssim h ^ {- \frac {1}{2}} \left(\| \sigma - \sigma_ {h} \| + \| \sigma - \tau_ {h} \| + \| u - u _ {h} \| + \| u - v _ {h} \|\right) + \\ h ^ {\frac {1}{2}} \left(\| \nabla_ {h} (\sigma - \tau_ {h}) \| + \| \nabla_ {h} (u - v _ {h}) \|\right), \\ \end{array}
$$

利用定理8.9和8.8可得

$$
\begin{array}{l} \| u - u _ {h} ^ {*} \| \lesssim h (\| q \| + \| \rho \|) + h \inf  _ {\tau_ {h} \in \Sigma_ {h}} \left(\| \sigma - \tau_ {h} \| + h \| \nabla_ {h} (\sigma - \tau_ {h}) \|\right) + \\ \inf  _ {v _ {h} \in V _ {h, 2}} \left(\| u - v _ {h} \| + h \| \nabla_ {h} (u - v _ {h}) \|\right). \\ \end{array}
$$

从而由引理8.10及有限元插值估计知本定理成立

# 8.3.5 变分形式II

本小节将HDG法改写为仅关于数值迹 $\hat{u}_h$ 的变分形式，从而数值实现时可以得到仅关于 $\hat{u}_h$ 的方程组，大大减少了自由度数.为了消掉 $u_{h}$ 和 $\sigma_h$ ，我们将(8.40)—(8.42）分解为两个单元 $K$ 上的子问题

第一个子问题为：任给 $m \in L^2(\mathcal{E}_h)$ ，求 $\mathcal{Q}_1 m \in V_h$ ， $\mathcal{Q}_2 m \in \Sigma_h$ 满足，对任意 $K \in \mathcal{M}_h$

$$
\begin{array}{l} \int_ {K} \kappa^ {- 1} Q _ {2} m \cdot \tau + \int_ {K} Q _ {1} m \nabla \cdot \tau = \int_ {\partial K} m \tau \cdot n _ {K}, \quad \forall \tau \in P _ {1} ^ {d} (K), (8.71) \\ \int_ {K} \mathcal {Q} _ {2} m \cdot \nabla v - \int_ {\partial K} \hat {\mathcal {Q}} m \cdot n _ {K} v = 0, \quad \forall v \in P _ {1} (K), (8.72) \\ \left. \hat {Q} m \right| _ {\partial K} = Q _ {2} m + \alpha (m - Q _ {1} m) n _ {K}. (8.73) \\ \end{array}
$$

第二个子问题为：任给 $f\in L^{2}(\varOmega)$ ，求 $\mathcal{Q}_1f\in V_h,\mathcal{Q}_2f\in \Sigma_h$ 满足 $\forall K\in \mathcal{M}_h$

$$
\begin{array}{l} \int_ {K} \kappa^ {- 1} Q _ {2} f \cdot \tau + \int_ {K} Q _ {1} f \nabla \cdot \tau = 0, \quad \forall \tau \in P _ {1} ^ {d} (K), (8.74) \\ \int_ {K} \mathcal {Q} _ {2} f \cdot \nabla v - \int_ {\partial K} \hat {\mathcal {Q}} f \cdot n _ {K} v = \int_ {K} f v, \quad \forall v \in P _ {1} (K), (8.75) \\ \left. \hat {\mathcal {Q}} f \right| _ {\partial K} = \mathcal {Q} _ {2} f - \alpha \mathcal {Q} _ {1} f n _ {K}. (8.76) \\ \end{array}
$$

显然，我们有

$$
u _ {h} = \mathcal {Q} _ {1} \hat {u} _ {h} + \mathcal {Q} _ {1} f; \quad \sigma_ {h} = \mathcal {Q} _ {2} \hat {u} _ {h} + \mathcal {Q} _ {2} f; \quad \hat {\sigma} _ {h} = \hat {\mathcal {Q}} \hat {u} _ {h} + \hat {\mathcal {Q}} f. \tag {8.77}
$$

由 (8.43) 知

$$
\langle \mu , [   [ \hat {\sigma} _ {h} ]   ] \rangle_ {\mathcal {E} _ {h} ^ {I N}} = 0, \quad \forall \mu \in S _ {h}.
$$

得 $\hat{u}_h$ 满足如下变分公式：求 $\hat{u}_h\in S_h$ 使得

$$
a _ {h} \left(\hat {u} _ {h}, \mu\right) := \left\langle \mu , \llbracket \hat {\mathcal {Q}} \hat {u} _ {h} \rrbracket \right\rangle_ {\mathcal {E} _ {h} ^ {I N}} = - \left\langle \mu , \llbracket \hat {\mathcal {Q}} f \rrbracket \right\rangle_ {\mathcal {E} _ {h} ^ {I N}}, \quad \forall \mu \in S _ {h}. \tag {8.78}
$$

# 8.3.6 习题

1. 给出 IIPG 法 $(\beta = 0)$ 的 $H^{1}$ 和 $L^{2}$ 误差估计  
2. 证明对任意 $v_{h} \in V_{h}$ , 存在 $v_{h}^{0} \in V_{h} \cap H^{1}(\Omega)$ , 满足

$$
\| v _ {h} - v _ {h} ^ {0} \| \lesssim \| h ^ {\frac {1}{2}} [ v _ {h} ] \| _ {\mathcal {E} _ {h} ^ {I}}, \quad \| \nabla (v _ {h} - v _ {h} ^ {0}) \| _ {\mathcal {M} _ {h}} \lesssim \| h ^ {- \frac {1}{2}} [ v _ {h} ] \| _ {\mathcal {E} _ {h} ^ {I}}.
$$

3. 假设网格 $\mathcal{M}_h$ 是协调的, $\gamma_1 = 0$ , 证明当 $\gamma_0 \to +\infty$ 时, SIPG 法的解收敛于有限元解 (提示: 利用题 2 结论).  
4. 对两点边值问题 $-u'' = 1, x \in (0,1), u(0) = u(1) = 0,$ 给出等距网格上的NIPG法（取 $\gamma_{1} = 0$ ），并作数值实验分别研究 $\gamma_0 = 1,0,h^{-2}$ 时的 $H^{1}$ 和 $L^2$ 误差阶  
5. 对LDG法，给出 $\| \kappa \nabla u - \sigma_h\|_{L^2 (\varOmega)}$ 的误差估计  
6. 证明 (8.78) 中的双线性形式满足：

$$
a _ {h} (\lambda , \mu) = \left(\kappa^ {- 1} \mathcal {Q} _ {2} \lambda , \mathcal {Q} _ {2} \mu\right) + \sum_ {K \in \mathcal {M} _ {h}} \alpha \left\langle \mathcal {Q} _ {1} \lambda - \lambda , \mathcal {Q} _ {1} \mu - \mu \right\rangle_ {\partial K}, \quad \forall \lambda , \mu \in S _ {h}.
$$

从而 $a_{h}$ 是对称的

![](images/0b49656fec8e9e8f6837e8f82fba5f0f2272a4f1ce3a19ebf97c80707ec092f6.jpg)

# 第9章

# 弱有限元法

弱有限元法[32, 33]是一种求解偏微分方程的高效数值方法. 弱有限元法的主要特点是引入弱函数作为近似函数, 并对弱函数定义弱微分算子. 弱函数包括两部分: 内部函数和边界函数. 内部函数通常采用间断的分片多项式, 边界函数是定义在单元边界上的多项式. 另外, 边界函数还用于单元与单元之间的联系. 因此, 我们在做离散的时候可以根据需要选择不同的内部函数空间和边界函数空间的组合. 对于这种完全间断的弱函数, 经典的微分算子不再适用. 因此, 我们根据Green公式来定义弱函数的弱微分算子, 进而替代变分形式中的经典微分算子. 另一特点是添加“稳定子”来保证内部函数与边界函数之间的弱连续性和稳定性.

弱有限元法的网格剖分是任意多边形或多面体剖分. 在接下来的分析中, 我们以三角形或四面体剖分为例, 并且考虑线性多项式的情形.

# 9.1 弱微分算子

# 9.1.1 广义弱微分算子

设 $\mathcal{T}_h$ 为有界区域 $\Omega \subset \mathbb{R}^d (d = 2,3)$ 的三角形或四面体剖分. 记 $\varepsilon_h$ 为剖分 $\mathcal{T}_h$ 中所有单元的边或面的集合. 对于任意的三角形或四面体单元 $K\in \mathcal{T}_h$ ，其边界记作 $\partial K$ . 令 $|K|$ 表示 $K$ 的面积, $h_K$ 表示单元 $K$ 的直径. $h = \max_{K\in \mathcal{T}_h}h_K$ 表示剖分 $\mathcal{T}_h$ 的网格尺寸. 单元 $K$ 上的弱函数（见图9.1）记作 $v = \{v_0,v_b\}$ ， $v_{0}\in L^{2}(K)$ ，且 $v_{b}\in L^{2}(\partial K)$ . 第一个分量 $v_{0}$ 表示 $v$ 在 $K$ 的内部的值，第二个分量 $v_{b}$ 表示 $v$ 在 $K$ 的边界上的值. 特别注意的是, $v_{b}$ 与 $v_{0}$ 在 $\partial K$ 上的迹没有必然联系. 记 $K$ 上所有弱函数构成的空间为 $W(K)$ ，即

$$
W (K) = \{v = \left\{v _ {0}, v _ {b} \right\}: v _ {0} \in L ^ {2} (K), v _ {b} \in L ^ {2} (\partial K) \}. \tag {9.1}
$$

![](images/e16ef5b7f80f4a3d3957bf48d56b9688d18c351059c6efcae4aa4836ba98aed9.jpg)  
图9.1 弱函数图像

这里我们强调 $v_{b}$ 在边 $e \in \mathcal{E}_{h}$ 上是单值函数。若 $v$ 在 $\Omega$ 上连续，则 $v = \{v, v\}$ 。记 $\pmb{n} = (n_{1}, n_{2}, \dots, n_{d})^{\mathrm{T}}$ 为边界 $\partial K$ 上的单位外法向量。记 $(\cdot, \cdot)_{K}$ 为空间 $L^{2}(K)$ 中的内积， $\langle \cdot, \cdot \rangle_{\partial K}$ 为空间 $L^{2}(\partial K)$ 中的内积。

定义9.1（一阶弱偏导数）对任意的 $v \in W(K)$ ， $v$ 的一阶弱偏导数 $\partial_{w,i} v$ ( $i = 1,2,\dots,d$ ) 定义为 $H^{1}(K)$ 上的有界线性泛函，满足

$$
\langle \partial_ {w, i} v, \phi \rangle_ {K} = - (v _ {0}, \partial_ {i} \phi) _ {K} + \langle v _ {b}, \phi n _ {i} \rangle_ {\partial K}, \quad \forall \phi \in H ^ {1} (K). \tag {9.2}
$$

定义9.2（弱梯度算子）对任意的 $v \in W(K)$ 中， $v$ 的弱梯度算子定义为向量值 Sobolev 空间 $[H^1(K)]^d$ 上的有界线性泛函，即

$$
\nabla_ {w} v = \left(\partial_ {w, 1} v, \partial_ {w, 2} v, \dots , \partial_ {w, d} v\right) ^ {\mathrm {T}}.
$$

那么，弱梯度 $\nabla_{w}v$ 满足

$$
\langle \nabla_ {w} v, q \rangle_ {K} = - (v _ {0}, \nabla \cdot q) _ {K} + \langle v _ {b}, q \cdot n \rangle_ {\partial K}, \quad \forall q \in [ H ^ {1} (K) ] ^ {d}. \tag {9.3}
$$

# 9.1.2 离散弱微分算子

对任意给定的非负整数 $r \geqslant 0$ , 记 $P_r(K)$ 为单元 $K$ 上次数不超过 $r$ 的多项式集合.

定义9.3（离散弱一阶偏导数）对任意的 $v \in W(K)$ ， $v$ 在 $x_{i}$ 方向的离散弱一阶偏导数 $\partial_{w,i,r} v$ ，定义为 $P_{r}(K)$ 中唯一的多项式，满足

$$
\left(\partial_ {w, i, r} v, \phi\right) _ {K} = - \left(v _ {0}, \partial_ {i} \phi\right) _ {K} + \left\langle v _ {b}, \phi n _ {i} \right\rangle_ {\partial K}, \quad \forall \phi \in P _ {r} (K). \tag {9.4}
$$

离散弱梯度算子定义为

$$
\nabla_ {w, r} v = \left(\partial_ {w, 1, r} v, \partial_ {w, 2, r} v, \dots , \partial_ {w, d, r} v\right) ^ {\mathrm {T}}. \tag {9.5}
$$

则满足

$$
\left(\nabla_ {w, r} v, \boldsymbol {q}\right) _ {K} = - \left(v _ {0}, \nabla \cdot \boldsymbol {q}\right) _ {K} + \left\langle v _ {b}, \boldsymbol {q} \cdot \boldsymbol {n} \right\rangle_ {\partial K}, \quad \forall \boldsymbol {q} \in \left[ P _ {r} (K) \right] ^ {d}. \tag {9.6}
$$

# 9.2 弱有限元数值格式

本节以二阶椭圆型方程为例详细介绍弱有限元数值格式

考虑 $\mathbb{R}^d (d = 2,3)$ 中有界区域 $\Omega$ 上的二阶椭圆型问题：

$$
- \nabla \cdot (\kappa \nabla u) = f, \quad \text {在} \Omega , \tag {9.7}
$$

$$
u = 0, \quad \text {在} \partial \Omega , \tag {9.8}
$$

其中右端函数 $f\in L^{2}(\varOmega)$ 给定， $\kappa = \kappa (x,y)$ 是区域 $\varOmega$ 上给定的对称正定矩阵值函数，即存在两个正常数 $\lambda_1,\lambda_2$ 使得

$$
\lambda_ {1} \boldsymbol {\xi} ^ {\mathrm {T}} \boldsymbol {\xi} \leqslant \boldsymbol {\xi} ^ {\mathrm {T}} \kappa \boldsymbol {\xi} \leqslant \lambda_ {2} \boldsymbol {\xi} ^ {\mathrm {T}} \boldsymbol {\xi}.
$$

此处, $\xi$ 是 $\mathbb{R}^d$ 空间中的列向量. 在接下来的部分, 我们考虑 $\kappa$ 是分片常数的情况.

对任意的 $v \in H_0^1(\Omega)$ , 通过对模型问题 (9.7) (9.8) 进行分部积分, 可以获得原始变量变分形式: 求 $u \in H_0^1(\Omega)$ , 使得

$$
(\kappa \nabla u, \nabla v) = (f, v), \quad \forall v \in H _ {0} ^ {1} (\Omega). \tag {9.9}
$$

在有限元剖分 $\mathcal{T}_h$ 上定义下列弱有限元空间 $V_{h}$ 为

$$
V _ {h} = \left\{v _ {h} = \left\{v _ {0}, v _ {b} \right\}: v _ {0} | _ {K} \in P _ {1} (K), v _ {b} | _ {e} \in P _ {0} (e), \forall e \subset \partial K, K \in \mathcal {T} _ {h} \right\}. \tag {9.10}
$$

那么，含有齐次边界条件的弱有限元空间可定义为

$$
V _ {h} ^ {0} = \left\{v _ {h} = \left\{v _ {0}, v _ {b} \right\}: v _ {h} \in V _ {h}, v _ {b} | _ {e} = 0, \quad \forall e \subset \partial \Omega \right\}. \tag {9.11}
$$

对于弱函数 $v_{h}$ ，简记 $\nabla_w v_h = \nabla_{w,0} v_h$ 为它的离散弱梯度(见(9.6)式).

对于 $K \in \mathcal{T}_h$ , 记 $Q_0$ 为从 $L^2(K)$ 到 $P_1(K)$ 的 $L^2$ 投影, $Q_b$ 为从 $L^2(e)$ 到 $P_0(e)$ 的 $L^2$ 投影. 记投影算子 $Q_h: H^1(\Omega) \to V_h$ , 使得在每个单元 $K \in \mathcal{T}_h$ 有

$$
Q _ {h} v = \left\{Q _ {0} v, Q _ {b} v \right\}. \tag {9.12}
$$

对于任意的 $w_{h},v_{h}\in V_{h}$ ，定义下列双线性形式：

$$
a \left(w _ {h}, v _ {h}\right) = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} w _ {h}, \nabla_ {w} v _ {h}\right) _ {K}, \tag {9.13}
$$

$$
s \left(w _ {h}, v _ {h}\right) = \sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \left\langle Q _ {b} w _ {0} - w _ {b}, Q _ {b} v _ {0} - v _ {b} \right\rangle_ {\partial K}, \tag {9.14}
$$

$$
a _ {s} \left(w _ {h}, v _ {h}\right) = a \left(w _ {h}, v _ {h}\right) + s \left(w _ {h}, v _ {h}\right). \tag {9.15}
$$

对于二阶椭圆型Dirichlet边值问题(9.7)(9.8)，相应的原始变量弱有限元格式为：求 $u_{h} = \{u_{0},u_{b}\} \in V_{h}^{0}$ ，使得

$$
a _ {s} \left(u _ {h}, v _ {h}\right) = (f, v _ {0}), \quad \forall v _ {h} = \left\{v _ {0}, v _ {b} \right\} \in V _ {h} ^ {0}. \tag {9.16}
$$

注9.1 在格式中添加稳定子 $s(u_h, v_h)$ 是为了保证数值解的弱连续性及稳定性.

# 9.2.1 适定性

为了证明弱有限元格式的适定性，我们在 $V_{h}$ 空间中定义如下半范数：

$$
\| v \| ^ {2} = a _ {s} (v, v), \quad \forall v \in V _ {h}. \tag {9.17}
$$

引理9.1 $\| \cdot \|$ 是 $V_h^0$ 空间上的范数

证明 令 $\|v\| = 0$ ，那么根据 $a_s(v, v)$ 的定义可以得到：

$$
\nabla_ {w} v = 0, \text {在} K \in \mathcal {T} _ {h};
$$

$$
Q _ {b} v _ {0} = v _ {b}, \text {在} e \subset \partial K.
$$

对任意的 $q \in [P_0(K)]^d$ ，根据离散弱梯度的定义 (9.6)，分部积分以及投影算子 $Q_b$ 的定义可以得到

$$
\begin{array}{l} 0 = (\nabla_ {w} v, q) _ {K} \\ = - \left(v _ {0}, \nabla \cdot \boldsymbol {q}\right) _ {K} + \left\langle v _ {b}, \boldsymbol {q} \cdot \boldsymbol {n} \right\rangle_ {\partial K} \\ = (\nabla v _ {0}, \boldsymbol {q}) _ {K} - \langle Q _ {b} v _ {0} - v _ {b}, \boldsymbol {q} \cdot \boldsymbol {n} \rangle_ {\partial K} \\ = (\nabla v _ {0}, \boldsymbol {q}) _ {K}. \\ \end{array}
$$

在上式中取 $q = \nabla v_{0}$ 可以得到，在每个单元 $K$ 上有 $\nabla v_{0} = 0$ 因此， $v_{0}$ 在每个单元 $K$ 上是常数. 又 $v_{0} = Q_{b}v_{0} = v_{b}$ ，故 $v_{0}$ 在整个区域 $\Omega$ 上为一个常数. 再根据 $v_{b}|_{\partial \Omega} = 0$ 故 $v_{0} = v_{b} = 0.$ 引理得证. □

引理9.2 对任意的 $v, w \in V_h$ ，有

$$
\begin{array}{l} \left| a _ {s} (v, w) \right| \leqslant \| | v | \| \| | w | \|, \\ a _ {s} (v, v) = \| | v \| ^ {2}. \\ \end{array}
$$

定理9.1 弱有限元数值格式（9.16）存在唯一解

证明 由于方程组（9.16）中的未知数个数与方程的个数相等，故解的存在性和唯一性等价。证明解的唯一性等价于证明齐次问题 $f = 0$ 时只有零解。当 $f = 0$ 时，取 $v_{h} = u_{h} \in V_{h}^{0}$ ，弱有限元格式（9.16）可写为

$$
a _ {s} \left(u _ {h}, u _ {h}\right) = \| \left| u _ {h} \right| \| ^ {2} = 0.
$$

根据引理9.1可知， $\| \cdot \|$ 是 $V_{h}^{0}$ 中的范数,故 $u_{h} = 0.$ 定理得证

# 9.2.2 $L^2$ 投影的误差估计

对 $K\in \mathcal{T}_h$ ，记 $\mathbb{Q}_h$ 是从 $[L^2 (K)]^d$ 到局部离散弱梯度空间 $[P_0(K)]^d$ 的 $L^2$ 投影.下面给出在弱有限元方法的误差估计中用到的一些重要不等式

引理9.3（投影不等式） 假设 $\mathcal{T}_h$ 是区域 $\varOmega$ 的形状正则剖分.那么，对于 $k = 0,1$ 及任意的 $\phi \in H^{k + 1}(\varOmega)$ 有

$$
\sum_ {K \in \mathcal {T} _ {h}} \| \phi - Q _ {0} \phi \| _ {K} ^ {2} + \sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {2} \| \nabla (\phi - Q _ {0} \phi) \| _ {K} ^ {2} \leqslant C h ^ {2 (k + 1)} \| \phi \| _ {k + 1} ^ {2}, \tag {9.18}
$$

$$
\sum_ {K \in \mathcal {T} _ {h}} \| \nabla \phi - \mathbb {Q} _ {h} (\nabla \phi) \| _ {K} ^ {2} \leqslant C h ^ {2} \| \phi \| _ {2} ^ {2}. \tag {9.19}
$$

其中， $C$ 是与网格尺寸 $h$ 和估计中的函数 $\phi$ 无关的正常数

# 9.2.3 $H^{1}$ 误差估计

令 $u_{h} = \{u_{0}, u_{b}\}$ 是弱有限元数值格式 (9.16) 的数值解, $u$ 是问题 (9.7) (9.8) 的精确解. 那么, 数值解 $u_{h}$ 与精确解 $u$ 的投影之间的误差记为

$$
e _ {h} = \left\{e _ {0}, e _ {b} \right\} = Q _ {h} u - u _ {h} = \left\{Q _ {0} u - u _ {0}, Q _ {b} u - u _ {b} \right\}.
$$

引理9.4 在单元 $K\in \mathcal{T}_h$ 上，有下列性质成立：

$$
\nabla_ {w} \left(Q _ {h} \phi\right) = \mathbb {Q} _ {h} (\nabla \phi), \quad \forall \phi \in H ^ {1} (K). \tag {9.20}
$$

证明 对任意的 $q \in [P_0(K)]^d$ ，根据离散弱梯度的定义 (9.6)，分部积分以及投影算子 $Q_h$ 和 $\mathbb{Q}_h$ 的定义可得，

$$
\begin{array}{l} \left(\nabla_ {w} \left(Q _ {h} \phi\right), q\right) _ {K} = - \left(Q _ {0} \phi , \nabla \cdot q\right) _ {K} + \left\langle Q _ {b} \phi , q \cdot n \right\rangle_ {\partial K} \\ = - (\phi , \nabla \cdot \boldsymbol {q}) _ {K} + \langle \phi , \boldsymbol {q} \cdot \boldsymbol {n} \rangle_ {\partial K} \\ = (\nabla \phi , \boldsymbol {q}) _ {K} = (\mathbb {Q} _ {h} (\nabla \phi), \boldsymbol {q}) _ {K}. \\ \end{array}
$$

在上式中取 $q = \nabla_{w}(Q_{h}\phi) - \mathbb{Q}_{h}(\nabla \phi)$ 便可得到(9.20).引理得证

引理9.5 对于 $v = \{v_0,v_b\} \in V_h^0$ ，我们有

$$
\sum_ {K \in \mathcal {T} _ {h}} \| \nabla v _ {0} \| _ {K} ^ {2} \leqslant C \| v \| ^ {2}. \tag {9.21}
$$

证明 对于 $v = \{v_0, v_b\} \in V_h^0$ ，根据离散弱梯度的定义（9.6）和投影算子 $Q_b$ 的定义可以得到

$$
\begin{array}{l} \left(\nabla_ {w} v, \boldsymbol {q}\right) _ {K} = - \left(v _ {0}, \nabla \cdot \boldsymbol {q}\right) _ {K} + \left\langle v _ {b}, \boldsymbol {q} \cdot \boldsymbol {n} \right\rangle_ {\partial K} \\ = (\nabla v _ {0}, \boldsymbol {q}) _ {K} + \langle v _ {b} - Q _ {b} v _ {0}, \boldsymbol {q} \cdot \boldsymbol {n} \rangle_ {\partial K}, \quad \forall \boldsymbol {q} \in [ P _ {0} (K) ] ^ {d}. \tag {9.22} \\ \end{array}
$$

在 (9.22) 中令 $q = \nabla v_0$ 可得

$$
\left(\nabla_ {w} v, \nabla v _ {0}\right) _ {K} = \left(\nabla v _ {0}, \nabla v _ {0}\right) _ {K} + \left\langle v _ {b} - Q _ {b} v _ {0}, \nabla v _ {0} \cdot \boldsymbol {n} \right\rangle_ {\partial K}.
$$

根据Cauchy-Schwarz不等式，迹不等式(6.77)和逆不等式(6.76）可以得到

$$
\| \nabla v _ {0} \| _ {K} ^ {2} \leqslant \| \nabla_ {w} v \| _ {K} \| \nabla v _ {0} \| _ {K} + C h _ {K} ^ {- \frac {1}{2}} \| v _ {b} - Q _ {b} v _ {0} \| _ {\partial K} \| \nabla v _ {0} \| _ {K}.
$$

因此, 我们有

$$
\left\| \nabla v _ {0} \right\| _ {K} \leqslant C \left(\left\| \nabla_ {w} v \right\| _ {K} ^ {2} + C h _ {K} ^ {- 1} \| v _ {b} - Q _ {b} v _ {0} \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}}.
$$

对单元 $K$ 进行求平方和可得(9.21)成立.引理得证

引理9.6 对任意的 $v \in V_h^0$ ，误差 $e_h$ 满足

$$
a _ {s} \left(e _ {h}, v\right) = \ell_ {u} (v) + s \left(Q _ {h} u, v\right), \tag {9.23}
$$

其中

$$
\ell_ {u} (v) = \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot \boldsymbol {n}, v _ {0} - v _ {b} \right\rangle_ {\partial K}. \tag {9.24}
$$

证明 在方程 (9.7) 的两端与 $v = \{v_0, v_b\}$ 中的 $v_0$ 作内积，并利用分部积分可得：

$$
\sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla u, \nabla v _ {0}) _ {K} - \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa \nabla u \cdot \boldsymbol {n}, v _ {0} - v _ {b} \right\rangle_ {\partial K} = (f, v _ {0}), \tag {9.25}
$$

这里, 我们用到了 $\sum_{K\in \mathcal{T}_h}\langle v_b,\pmb {\kappa}\nabla u\cdot \pmb {n}\rangle_{\partial K} = 0.$ 对于 $u\in H^{2}(\varOmega)$ 和 $v\in V_h^0$ ，根据(9.20)，离散弱梯度的定义(9.6)，分部积分和投影算子 $\mathbb{Q}_h$ 的定义可得：

$$
\begin{array}{l} \left(\kappa \nabla_ {w} Q _ {h} u, \nabla_ {w} v\right) _ {K} = \left(\kappa \mathbb {Q} _ {h} (\nabla u), \nabla_ {w} v\right) _ {K} \\ = - \left(v _ {0}, \nabla \cdot (\kappa \mathbb {Q} _ {h} (\nabla u))\right) _ {K} + \left\langle v _ {b}, (\kappa \mathbb {Q} _ {h} (\nabla u)) \cdot \boldsymbol {n} \right\rangle_ {\partial K} \\ = (\nabla v _ {0}, \kappa \mathbb {Q} _ {h} (\nabla u)) _ {K} - \langle v _ {0} - v _ {b}, (\kappa \mathbb {Q} _ {h} (\nabla u)) \cdot \boldsymbol {n} \rangle_ {\partial K} \\ = (\nabla v _ {0}, \kappa \nabla u) _ {K} - \langle v _ {0} - v _ {b}, (\kappa \mathbb {Q} _ {h} (\nabla u)) \cdot \boldsymbol {n} \rangle_ {\partial K}. \tag {9.26} \\ \end{array}
$$

将 (9.25) 代入 (9.26) 中可得

$$
\begin{array}{l} \sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla_ {w} Q _ {h} u, \nabla_ {w} v) _ {K} = (f, v _ {0}) + \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot n, v _ {0} - v _ {b} \right\rangle_ {\partial K} \\ = (f, v _ {0}) + \bar {\ell_ {u}} (v). \\ \end{array}
$$

在上式两端添加稳定子 $s(Q_h u, v)$ 可得

$$
a _ {s} \left(Q _ {h} u, v\right) = \left(f, v _ {0}\right) + \ell_ {u} (v) + s \left(Q _ {h} u, v\right). \tag {9.27}
$$

(9.27) 减去 (9.16) 可得

$$
a _ {s} \left(e _ {h}, v\right) = \ell_ {u} (v) + s \left(Q _ {h} u, v\right).
$$

引理得证.

引理9.7 对任意的 $w \in H^{2}(\Omega)$ 和 $v = \{v_{0}, v_{b}\} \in V_{h}$ , 我们有

$$
\left| \ell_ {w} (v) \right| \leqslant C h \| w \| _ {2} \| v \|, \tag {9.28}
$$

$$
| s (Q _ {h} w, v) | \leqslant C h \| w \| _ {2} \| v \| \|. \tag {9.29}
$$

证明 根据Cauchy-Schwarz不等式，迹不等式(6.77)和投影不等式(9.19)可得：

$$
\begin{array}{l} \left| \ell_ {w} (v) \right| = \left| \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla w - \mathbb {Q} _ {h} \nabla w) \cdot \boldsymbol {n}, v _ {0} - v _ {b} \right\rangle_ {\partial K} \right| \\ \leqslant C \sum_ {K \in \mathcal {T} _ {h}} \| \kappa (\nabla w - \mathbb {Q} _ {h} \nabla w) \| _ {\partial K} \| v _ {0} - v _ {b} \| _ {\partial K} \\ \leqslant C \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} \| \kappa (\nabla w - \mathbb {Q} _ {h} \nabla w) \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| v _ {0} - v _ {b} \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \\ \leqslant C h \| w \| _ {2} \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| v _ {0} - v _ {b} \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}}. \tag {9.30} \\ \end{array}
$$

根据迹不等式 (6.77) 和投影不等式 (9.18) 可以得到

$$
\begin{array}{l} \| v _ {0} - v _ {b} \| _ {\partial K} \leqslant \| v _ {0} - Q _ {b} v _ {0} \| _ {\partial K} + \| Q _ {b} v _ {0} - v _ {b} \| _ {\partial K} \\ \leqslant C h _ {K} ^ {\frac {1}{2}} \| \nabla v _ {0} \| _ {K} + \| Q _ {b} v _ {0} - v _ {b} \| _ {\partial K}. \\ \end{array}
$$

将上式代入估计 (9.30), 由引理 9.5 可以得到

$$
\left| \ell_ {w} (v) \right| \leqslant C h \| w \| _ {2} \left(\sum_ {K \in \mathcal {T} _ {h}} \left(\| \nabla v _ {0} \| _ {K} ^ {2} + h _ {K} ^ {- 1} \| Q _ {b} v _ {0} - v _ {b} \| _ {\partial K} ^ {2}\right)\right) ^ {\frac {1}{2}} \leqslant C h \| w \| _ {2} \| v \| \|.
$$

关于估计 (9.29), 根据投影算子 $Q_{b}$ 的定义, Cauchy-Schwarz 不等式, 迹不等式 (6.77) 和投影不等式 (9.18) 可得

$$
\begin{array}{l} | s (Q _ {h} w, v) | = \left| \sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \langle Q _ {b} (Q _ {0} w) - Q _ {b} w, Q _ {b} v _ {0} - v _ {b} \rangle_ {\partial K} \right| \\ = \left| \sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \langle Q _ {0} w - w, Q _ {b} v _ {0} - v _ {b} \rangle_ {\partial K} \right| \\ \leqslant C \left(\sum_ {K \in \mathcal {T} _ {h}} \left(h _ {K} ^ {- 2} \| Q _ {0} w - w \| _ {K} ^ {2} + \| \nabla (Q _ {0} w - w) \| _ {K} ^ {2}\right)\right) ^ {\frac {1}{2}} \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| Q _ {b} v _ {0} - v _ {b} \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \\ \leqslant C h \| w \| _ {2} \| v \|. \\ \end{array}
$$

引理得证.

定理9.2 假设精确解 $u\in H^{2}(\varOmega)$ ，那么存在正常数 $C$ 使得

$$
\left\| u _ {h} - Q _ {h} u \right\| \leqslant C h \| u \| _ {2}. \tag {9.31}
$$

证明 在误差方程 (9.23) 中令 $v = e_h$ 可以得到

$$
\left\| \left| e _ {h} \right| \right\| ^ {2} = \ell_ {u} \left(e _ {h}\right) + s \left(Q _ {h} u, e _ {h}\right). \tag {9.32}
$$

再根据引理9.7可以得到

$$
\| \left\| e _ {h} \right\| \left\| \leqslant C h \right\| u \| _ {2}.
$$

定理得证

# 9.2.4 $L^2$ 误差估计

本节给出弱有限元解 $u_{h}$ 和精确解 $\mathcal{U}$ 之间的 $L^2$ 误差估计.我们主要通过对偶技巧来估计 $L^2$ 误差.考虑下列对偶问题：求 $\varPhi\in H_0^1(\varOmega)\cap H^2(\varOmega)$ 满足

$$
- \nabla \cdot (\kappa \nabla \Phi) = e _ {0} \text {, 在} \Omega . \tag {9.33}
$$

假设上述对偶问题有 $H^2$ 正则性，即存在正常数 $C$ 使得

$$
\| \Phi \| _ {2} \leqslant C \| e _ {0} \|. \tag {9.34}
$$

定理9.3 假设 $u_{h}\in V_{h}$ 是弱有限元数值解，精确解 $u\in H^{2}(\varOmega)$ ，那么存在一个正常数 $C$ 使得

$$
\| u - u _ {0} \| \leqslant C h ^ {2} \| u \| _ {2}. \tag {9.35}
$$

证明 在方程 (9.33) 两端与 $e_0$ 作内积, 并利用分部积分可得

$$
\begin{array}{l} \left\| e _ {0} \right\| ^ {2} = - (\nabla \cdot (\kappa \nabla \Phi), e _ {0}) \\ = \sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla \Phi , \nabla e _ {0}) _ {K} - \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa \nabla \Phi \cdot \boldsymbol {n}, e _ {0} - e _ {b} \right\rangle_ {\partial K}, \tag {9.36} \\ \end{array}
$$

这里我们用到了 $\sum_{K\in \mathcal{T}_h}\langle \kappa \nabla \Phi \cdot \pmb {n},e_b\rangle_{\partial K} = 0.$ 在（9.26）中，令 $u = \varPhi$ 和 $v = e_h$ 得到

$$
\left(\kappa \nabla_ {w} Q _ {h} \Phi , \nabla_ {w} e _ {h}\right) _ {K} = \left(\kappa \nabla \Phi , \nabla e _ {0}\right) _ {K} - \left\langle \kappa \mathbb {Q} _ {h} (\nabla \Phi) \cdot \boldsymbol {n}, e _ {0} - e _ {b} \right\rangle_ {\partial K}. \tag {9.37}
$$

将 (9.37) 代入 (9.36) 可以得到

$$
\begin{array}{l} \left\| e _ {0} \right\| ^ {2} = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} Q _ {h} \Phi , \nabla_ {w} e _ {h}\right) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa \left(\mathbb {Q} _ {h} (\nabla \Phi) - \nabla \Phi\right) \cdot \boldsymbol {n}, e _ {0} - e _ {b} \right\rangle_ {\partial K} \\ = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} Q _ {h} \Phi , \nabla_ {w} e _ {h}\right) _ {K} - \ell_ {\Phi} \left(e _ {h}\right). \tag {9.38} \\ \end{array}
$$

根据误差方程 (9.23) 可知

$$
\sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla_ {w} e _ {h}, \nabla_ {w} Q _ {h} \Phi) _ {K} = \ell_ {u} (Q _ {h} \Phi) + s (Q _ {h} u, Q _ {h} \Phi) - s (e _ {h}, Q _ {h} \Phi). \tag {9.39}
$$

由（9.38）（9.39）可得

$$
\left\| e _ {0} \right\| ^ {2} = \ell_ {u} \left(Q _ {h} \Phi\right) + s \left(Q _ {h} u, Q _ {h} \Phi\right) - s \left(e _ {h}, Q _ {h} \Phi\right) - \ell_ {\Phi} \left(e _ {h}\right). \tag {9.40}
$$

接下来估计 (9.40) 的每一项. 根据三角不等式可得

$$
\left| \ell_ {u} \left(Q _ {h} \Phi\right) \right| = \left| \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot n, Q _ {0} \Phi - Q _ {b} \Phi \right\rangle_ {\partial K} \right|
$$

$$
\begin{array}{l} \leqslant \left| \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot \boldsymbol {n}, Q _ {0} \Phi - \Phi \right\rangle_ {\partial K} \right| + \\ \left| \sum_ {K \in \mathcal {T} _ {h}} \langle \boldsymbol {\kappa} (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot \boldsymbol {n}, \Phi - Q _ {b} \Phi \rangle_ {\partial K} \right|. \tag {9.41} \\ \end{array}
$$

使用投影算子 $Q_{b}$ 的定义和在 $\partial \Omega$ 上 $\varPhi=0$ 可以得到

$$
\sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot \boldsymbol {n}, \Phi - Q _ {b} \Phi \right\rangle_ {\partial K} = \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa \nabla u \cdot \boldsymbol {n}, \Phi - Q _ {b} \Phi \right\rangle_ {\partial K} = 0. \tag {9.42}
$$

根据迹不等式 (6.77) 和投影不等式 (9.18) (9.19) 可得

$$
\left(\sum_ {K \in \mathcal {T} _ {h}} \| Q _ {0} \Phi - \Phi \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \leqslant C h ^ {\frac {3}{2}} \| \Phi \| _ {2}, \tag {9.43}
$$

和

$$
\left(\sum_ {K \in \mathcal {T} _ {h}} \| \kappa (\nabla u - \mathbb {Q} _ {h} \nabla u) \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \leqslant C h ^ {\frac {1}{2}} \| u \| _ {2}. \tag {9.44}
$$

根据Cauchy-Schwarz不等式和上述两个估计可以得到

$$
\begin{array}{l} \left| \sum_ {K \in \mathcal {T} _ {h}} \langle \boldsymbol {\kappa} (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot \boldsymbol {n}, Q _ {0} \Phi - \Phi \rangle_ {\partial K} \right| \\ \leqslant C \left(\sum_ {K \in \mathcal {T} _ {h}} \| \kappa (\nabla u - \mathbb {Q} _ {h} \nabla u) \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {K \in \mathcal {T} _ {h}} \| Q _ {0} \Phi - \Phi \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \\ \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \tag {9.45} \\ \end{array}
$$

再结合（9.42）和（9.45）可以得到

$$
\left| \ell_ {u} \left(Q _ {h} \Phi\right) \right| \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \tag {9.46}
$$

类似地, 根据投影算子 $Q_{b}$ 的定义, 迹不等式 (6.77) 和投影不等式 (9.18) 可得

$$
\begin{array}{l} \left| s \left(Q _ {h} u, Q _ {h} \Phi\right) \right| \leqslant \left| \sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \left\langle Q _ {b} \left(Q _ {0} u\right) - Q _ {b} u, Q _ {b} \left(Q _ {0} \Phi\right) - Q _ {b} \Phi \right\rangle_ {\partial K} \right| \\ \leqslant \sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| Q _ {b} (Q _ {0} u - u) \| _ {\partial K} \| Q _ {b} (Q _ {0} \Phi - \Phi) \| _ {\partial K} \\ \leqslant \sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| Q _ {0} u - u \| _ {\partial K} \| Q _ {0} \Phi - \Phi \| _ {\partial K} \\ \leqslant C \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| Q _ {0} u - u \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| Q _ {0} \Phi - \Phi \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \\ \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \tag {9.47} \\ \end{array}
$$

根据估计 (9.29) 和 (9.31) 可得

$$
\left| s \left(e _ {h}, Q _ {h} \Phi\right) \right| \leqslant C h \| \Phi \| _ {2} \| e _ {h} \| \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \tag {9.48}
$$

类似地，根据估计（9.28）和（9.31）可得

$$
\left| \ell_ {\Phi} \left(e _ {h}\right) \right| \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \tag {9.49}
$$

现在将估计(9.46)一(9.49)代入(9.40)可得

$$
\left\| e _ {0} \right\| ^ {2} \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \tag {9.50}
$$

再结合 $H^2$ 正则性估计 (9.34) 和三角不等式可以得到估计 (9.35). 定理得证

# 9.3 无稳定子弱有限元数值格式

本节给出求解二阶椭圆型问题 (9.7) (9.8) 的无稳定子弱有限元数值格式. 无稳定子弱有限元方法 [34], 通过提高弱梯度算子值域空间中多项式的次数, 去掉了数值格式中的稳定子, 简化了数值格式.

在有限元剖分 $\mathcal{T}_h$ 上定义下列弱有限元空间 $V_h$

$$
V _ {h} = \left\{v _ {h} = \left\{v _ {0}, v _ {b} \right\}: v _ {0} | _ {K} \in P _ {1} (K), v _ {b} | _ {e} \in P _ {1} (e), \forall e \subset \partial K, K \in \mathcal {T} _ {h} \right\}.
$$

那么，含有齐次边界条件的弱有限元空间可定义为

$$
V _ {h} ^ {0} = \left\{v _ {h}: v _ {h} \in V _ {h}, v _ {b} | _ {e} = 0, e \subset \partial \Omega \right\}.
$$

在无稳定子弱有限元法中, 我们将弱函数 $v$ 的离散弱梯度的定义延伸至 $V_{h} + H^{1}(\varOmega)$

定义9.4对任意的 $v\in V_h + H^1 (\Omega)$ ，在每个单元 $K\in \mathcal{T}_h$ 定义离散弱梯度 $(\nabla_w v)|_K\in [P_j(K)]^d (j > 1)$ ，满足

$$
\left(\nabla_ {w} v, \boldsymbol {q}\right) _ {K} = - \left(v _ {0}, \nabla \cdot \boldsymbol {q}\right) _ {K} + \left\langle v _ {b}, \boldsymbol {q} \cdot \boldsymbol {n} \right\rangle_ {\partial K}, \quad \forall \boldsymbol {q} \in [ P _ {j} (K) ] ^ {d}. \tag {9.51}
$$

在下一小节的证明中可以看到，弱梯度的次数可以选取为 $j = n$ ，这里 $n$ 表示多边形单元的边数或多面体单元的面数。

对于 $K \in \mathcal{T}_h$ , 记 $Q_0$ 为从 $L^2(K)$ 到 $P_1(K)$ 的 $L^2$ 投影, $Q_b$ 为从 $L^2(e)$ 到 $P_1(e)$ 的 $L^2$ 投影. 记投影算子 $Q_h: H^1(\Omega) \to V_h$ , 使得在每个单元 $K \in \mathcal{T}_h$ , 有

$$
Q _ {h} v = \left\{Q _ {0} v, Q _ {b} v \right\}. \tag {9.52}
$$

因此，对于二阶椭圆型Dirichlet边值问题(9.7)(9.8)，相应的无稳定子弱有限元数值格式为：求 $u_{h} = \{u_{0},u_{b}\} \in V_{h}^{0}$ 满足方程

$$
\sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla_ {w} u _ {h}, \nabla_ {w} v) _ {K} = \sum_ {K \in \mathcal {T} _ {h}} (f, v _ {0}) _ {K}, \quad \forall v _ {h} = \{v _ {0}, v _ {b} \} \in V _ {h} ^ {0}. \tag {9.53}
$$

# 9.3.1 适定性

对于任意的 $v \in V_h + H^1(\Omega)$ , 定义下列半范数:

$$
\| v \| ^ {2} = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} v, \nabla_ {w} v\right) _ {K}. \tag {9.54}
$$

引理9.8 设 $K$ 是一个直径为 $h_K$ 的凸多边形或凸多面体，边数或面数为 $n$ ，记 $e, e_2, \dots, e_n$ 为它的边。对于任意给定的多项式函数 $q_0 \in P_1(e)$ ，定义多项式 $q = \lambda_2\lambda_3 \dots \lambda_n q^* \in P_n(K)$ ，且满足

$$
\langle q - q _ {0}, p \rangle_ {e} = 0, \quad \forall p \in P _ {1} (e), \tag {9.55}
$$

$$
(q, r) _ {K} = 0, \quad \forall r \in P _ {0} (K), \tag {9.56}
$$

其中， $\lambda_{i}\in P_{1}(K)$ 在 $e_i$ $(i = 2,3,\dots ,n)$ 上为0，并假设在 $\mathcal{E}$ 的重心处值为 $1,q^{*}\in P_{1}(K)$ 则存在与单元 $K$ 和 $q_{0}$ 无关的常数 $C$ ，使得以下估计成立

$$
\left\| q \right\| _ {K} \leqslant C h _ {\hat {K}} ^ {\frac {1}{2}} \left\| q _ {0} \right\| _ {e}.
$$

引理9.9 对任意的 $v = \{v_0,v_b\} \in V_h$ 和单元 $K\in \mathcal{T}_h$ ，存在正常数 $C$ 使得如下不等式成立

$$
\left\| v _ {b} - v _ {0} \right\| _ {\partial K} ^ {2} \leqslant C h _ {K} \left\| \nabla_ {w} v \right\| _ {K} ^ {2}. \tag {9.57}
$$

证明 对任意 $v = \{v_0, v_b\} \in V_h$ ，根据弱梯度的定义（9.51）和分部积分，可以得到

$$
\left(\nabla_ {w} v, \boldsymbol {q}\right) _ {K} = \left(\nabla v _ {0}, \boldsymbol {q}\right) _ {K} + \left\langle v _ {b} - v _ {0}, \boldsymbol {q} \cdot \boldsymbol {n} \right\rangle_ {\partial K}, \quad \forall \boldsymbol {q} \in [ P _ {j} (K) ] ^ {d}. \tag {9.58}
$$

由引理9.8知存在 $q_0\in [P_n(K)]^d$ ，满足

$$
\left(\nabla v _ {0}, \boldsymbol {q} _ {0}\right) _ {K} = 0, \quad \langle v _ {b} - v _ {0}, \boldsymbol {q} _ {0} \cdot \boldsymbol {n} \rangle_ {\partial K \backslash e} = 0, \quad \langle v _ {b} - v _ {0}, \boldsymbol {q} _ {0} \cdot \boldsymbol {n} \rangle_ {e} = \| v _ {0} - v _ {b} \| _ {e} ^ {2},
$$

和

$$
\left\| \boldsymbol {q} _ {0} \right\| _ {K} \leqslant C h _ {K} ^ {\frac {1}{2}} \left\| v _ {b} - v _ {0} \right\| _ {e}. \tag {9.59}
$$

在 (9.58) 取 $\pmb{q} = \pmb{q}_0$ ，有

$$
\left(\nabla_ {w} v, \boldsymbol {q} _ {0}\right) _ {K} = \left\| v _ {b} - v _ {0} \right\| _ {e} ^ {2}.
$$

由Cauchy-Schwarz不等式和（9.59）可得

$$
\left\| v _ {b} - v _ {0} \right\| _ {e} ^ {2} \leqslant C \left\| \nabla_ {w} v \right\| _ {K} \left\| q _ {0} \right\| _ {K} \leqslant C h _ {K} ^ {\frac {1}{2}} \left\| \nabla_ {w} v \right\| _ {K} \left\| v _ {0} - v _ {b} \right\| _ {e},
$$

这表明

$$
h _ {K} ^ {- \frac {1}{2}} \| v _ {0} - v _ {b} \| _ {\partial K} \leqslant C \| \nabla_ {w} v \| _ {K}. \tag {9.60}
$$

引理得证.

从上述引理的证明中可以看出，为了确保估计（9.57）成立，我们选取弱梯度的次数为 $j = n = d + 1$

引理9.10 $\| \cdot \|$ 是 $V_{h}^{0}$ 空间上的范数

证明 假设 $v \in V_h^0$ 且 $\| v \| = 0$ . 根据 $\|\cdot\|$ 的定义可得

$$
\| v \| ^ {2} = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} v, \nabla_ {w} v\right) _ {K} = 0,
$$

即

$$
\left. \nabla_ {w} v \right| _ {K} = 0. \tag {9.61}
$$

根据估计 (9.57) 可得

$$
\left\| v _ {0} - v _ {b} \right\| _ {\partial K} = 0. \tag {9.62}
$$

对任意的单元 $K \in \mathcal{T}_h$ 和 $\pmb{q} \in [P_n(K)]^d$ , 根据 (9.58) 和 (9.61) (9.62) 可得

$$
(\nabla v _ {0}, \boldsymbol {q}) _ {K} = (\nabla_ {w} v, \boldsymbol {q}) _ {K} + \langle v _ {0} - v _ {b}, \boldsymbol {q} \cdot \boldsymbol {n} \rangle_ {\partial K} = 0.
$$

取 $q = \nabla v_{0}$ 可得

$$
\| \nabla v _ {0} \| _ {K} = 0.
$$

这意味着 $v_{0}$ 在每个单元 $K \in \mathcal{T}_h$ 上是一个常数. 又因为 $v_{b}$ 在边界 $\partial \Omega$ 上为 0. 再结合 (9.62) 可得, 在整个区域 $\Omega$ 上有 $v_{0} = v_{b} = 0$ . 引理得证. □

引理9.11 无稳定子弱有限元数值格式（9.53）存在唯一解

# 9.3.2 $H^{1}$ 误差估计

令 $u_{h} = \{u_{0}, u_{b}\}$ 是弱有限元数值格式 (9.53) 的数值解, $u$ 为问题 (9.7) (9.8) 的精确解. 那么数值解 $u_{h}$ 与精确解 $u$ 之间的误差 $e_{h}$ 记为

$$
e _ {h} = \left\{e _ {0}, e _ {b} \right\} = u - u _ {h} = \left\{u - u _ {0}, u - u _ {b} \right\}.
$$

另外，记 $\varepsilon_{h} = Q_{h}u - u_{h}$ ，那么误差 $e_h$ 可记为

$$
e _ {h} = (u - Q _ {h} u) + \varepsilon_ {h}.
$$

引理9.12 在单元 $K \in \mathcal{T}_h$ 上，记 $\mathbb{Q}_h$ 是到局部离散弱梯度空间 $[P_n(K)]^d$ 的 $L^2$

投影，则有下列性质成立

$$
\nabla_ {w} \phi = \mathbb {Q} _ {h} \nabla \phi , \quad \forall \phi \in H ^ {1} (K). \tag {9.63}
$$

证明 对任意 $q \in [P_n(K)]^d$ , 根据 (9.58), 可以得到

$$
\left(\nabla_ {w} \phi , \boldsymbol {q}\right) _ {K} = \left(\nabla \phi , \boldsymbol {q}\right) _ {K} = \left(\mathbb {Q} _ {h} \nabla \phi , \boldsymbol {q}\right) _ {K},
$$

在上式中取 $q = \nabla_{w}\phi -\mathbb{Q}_{h}\nabla \phi$ 可得式（9.63）成立

引理9.13 对任意的 $v \in V_h^0$ ，误差 $e_h$ 满足

$$
\sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla_ {w} e _ {h}, \nabla_ {w} v) _ {K} = \ell_ {u} (v), \tag {9.64}
$$

其中

$$
\ell_ {u} (v) = \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot \boldsymbol {n}, v _ {0} - v _ {b} \right\rangle_ {\partial K}.
$$

证明 在方程 (9.7) 的两端与 $v = \{v_0, v_b\}$ 中的 $v_0$ 作内积，并利用分部积分可得

$$
\sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla u, \nabla v _ {0}) _ {K} - \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa \nabla u \cdot n, v _ {0} - v _ {b} \right\rangle_ {\partial K} = \sum_ {K \in \mathcal {T} _ {h}} (f, v _ {0}) _ {K}. \tag {9.65}
$$

这里我们用到 $\sum_{K\in \mathcal{T}_h}\langle \kappa \nabla u\cdot n,v_b\rangle_{\partial K} = 0$ 根据(9.63)，分部积分和弱梯度的定义(9.51)可得

$$
\begin{array}{l} \sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla u, \nabla v _ {0}) _ {K} = \sum_ {K \in \mathcal {T} _ {h}} (\kappa \mathbb {Q} _ {h} \nabla u, \nabla v _ {0}) _ {K} \\ = \sum_ {K \in \mathcal {T} _ {h}} - (v _ {0}, \nabla \cdot (\kappa \mathbb {Q} _ {h} \nabla u)) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \langle v _ {0}, \kappa \mathbb {Q} _ {h} \nabla u \cdot \boldsymbol {n} \rangle_ {\partial K} \\ = \sum_ {K \in \mathcal {T} _ {h}} (\kappa \mathbb {Q} _ {h} \nabla u, \nabla_ {w} v) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \left\langle v _ {0} - v _ {b}, \kappa \mathbb {Q} _ {h} \nabla u \cdot n \right\rangle_ {\partial K} \\ = \sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla_ {w} u, \nabla_ {w} v) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \langle v _ {0} - v _ {b}, \kappa \mathbb {Q} _ {h} \nabla u \cdot n \rangle_ {\partial K} \tag {9.66} \\ \end{array}
$$

结合（9.65）和（9.66）可以得到

$$
\sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla_ {w} u, \nabla_ {w} v) _ {K} = \sum_ {K \in \mathcal {T} _ {h}} (f, v _ {0}) _ {K} + \ell_ {u} (v). \tag {9.67}
$$

(9.67) 减去数值格式 (9.53) 得到

$$
\sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla_ {w} e _ {h}, \nabla_ {w} v) _ {K} = \ell_ {u} (v), \quad \forall v \in V _ {h} ^ {0}.
$$

引理得证

引理9.14 对任意的 $w \in H^{2}(\Omega)$ 和 $v = \{v_{0}, v_{b}\} \in V_{h}^{0}$ , 有下列估计成立:

$$
\left| \ell_ {w} (v) \right| \leqslant C h \| w \| _ {2} \| v \| \tag {9.68}
$$

证明 根据 Cauchy-Schwarz 不等式, 迹不等式 (6.77), 投影不等式 (9.19) 和估计 (9.57), 我们有

$$
\begin{array}{l} | \ell_ {w} (v) | = \left| \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla w - \mathbb {Q} _ {h} \nabla w) \cdot \boldsymbol {n}, v _ {0} - v _ {b} \right\rangle_ {\partial K} \right| \\ \leqslant C \sum_ {K \in \mathcal {T} _ {h}} \| \nabla w - \mathbb {Q} _ {h} \nabla w \| _ {\partial K} \| v _ {0} - v _ {b} \| _ {\partial K} \\ \leqslant C \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} \| \nabla w - \mathbb {Q} _ {h} \nabla w \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| v _ {0} - v _ {b} \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \\ \leqslant C h \| w \| _ {2} \| v \|. \\ \end{array}
$$

引理得证.

引理9.15 对任意的 $w \in H^2(\Omega)$ , 存在正常数 $C$ 使得下列估计成立

$$
\| \left| w - Q _ {h} w \right| \| \leqslant C h \| w \| _ {2}. \tag {9.69}
$$

证明 对任意的 $q \in [P_n(K)]^d$ ，根据弱梯度的定义 (9.51)，分部积分，迹不等式 (6.77) 和逆不等式 (6.76) 可得

$$
\begin{array}{l} \left(\nabla_ {w} (w - Q _ {h} w), \boldsymbol {q}\right) _ {K} = - \left(w - Q _ {0} w, \nabla \cdot \boldsymbol {q}\right) _ {K} + \langle w - Q _ {b} w, \boldsymbol {q} \cdot \boldsymbol {n} \rangle_ {\partial K} \\ = (\nabla (w - Q _ {0} w), \boldsymbol {q}) _ {K} + \langle Q _ {0} w - Q _ {b} w, \boldsymbol {q} \cdot \boldsymbol {n} \rangle_ {\partial K} \\ \leqslant C \| \nabla (w - Q _ {0} w) \| _ {K} \| q \| _ {K} + C h _ {K} ^ {- \frac {1}{2}} \| w - Q _ {0} w \| _ {\partial K} \| q \| _ {K} \\ \leqslant C h \| w \| _ {2, K} \| \boldsymbol {q} \| _ {K}. \\ \end{array}
$$

在上式中令 $q = \nabla_w(w - Q_h w)$ 并对 $K$ 求和可得

$$
\| \left| w - Q _ {h} w \right| \| \leqslant C h \| w \| _ {2}.
$$

引理得证.

定理9.4 假设精确解 $u \in H^2(\Omega)$ ，那么存在正常数 $C$ 使得

$$
\| u - u _ {h} \| \leqslant C h \| u \| _ {2}. \tag {9.70}
$$

证明 根据 $\| \cdot \|$ 的定义可得

$$
\begin{array}{l} \| \left| e _ {h} \right| \| ^ {2} = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} e _ {h}, \nabla_ {w} e _ {h}\right) _ {K} \\ = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \left(\nabla_ {w} u - \nabla_ {w} u _ {h}\right), \nabla_ {w} e _ {h}\right) _ {K} \\ \end{array}
$$

$$
\begin{array}{l} = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \left(\nabla_ {w} Q _ {h} u - \nabla_ {w} u _ {h}\right), \nabla_ {w} e _ {h}\right) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \left(\nabla_ {w} u - \nabla_ {w} Q _ {h} u\right), \nabla_ {w} e _ {h}\right) _ {K} \\ = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} e _ {h}, \nabla_ {w} \varepsilon_ {h}\right) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \left(\nabla_ {w} u - \nabla_ {w} Q _ {h} u\right), \nabla_ {w} e _ {h}\right) _ {K}. \tag {9.71} \\ \end{array}
$$

接下来我们估计 (9.71) 的每一项．在误差方程 (9.64) 中令 $v = \varepsilon_h \in V_h^0$ 并根据估计 (9.68)一(9.69)和Young不等式可得

$$
\begin{array}{l} \left| \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} e _ {h}, \nabla_ {w} \varepsilon_ {h}\right) _ {K} \right| = | \ell_ {u} (\varepsilon_ {h}) | \\ \leqslant C h \| u \| _ {2} \| \varepsilon_ {h} \| \\ \leqslant C h \| u \| _ {2} \| Q _ {h} u - u _ {h} \| \\ \leqslant C h \| u \| _ {2} \left(\| Q _ {h} u - u \| + \| u - u _ {h} \|\right) \\ \leqslant C h ^ {2} \| u \| _ {2} ^ {2} + \frac {1}{4} \| e _ {h} \| ^ {2}. \tag {9.72} \\ \end{array}
$$

类似地，根据估计（9.69）和Young不等式可得

$$
\begin{array}{l} \left| \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \left(\nabla_ {w} u - \nabla_ {w} Q _ {h} u\right), \nabla_ {w} e _ {h}\right) _ {K} \right| \leqslant C \| | u - Q _ {h} u | \| \| e _ {h} \| \\ \leqslant C h ^ {2} \| u \| _ {2} ^ {2} + \frac {1}{4} \| e _ {h} \| ^ {2}. \tag {9.73} \\ \end{array}
$$

结合估计 (9.71)一(9.73)可得

$$
\| e _ {h} \| \leqslant C h \| u \| _ {2}.
$$

定理得证

# 9.3.3 $L^2$ 误差估计

接下来我们使用对偶方法来得到 $L^2$ 误差估计. 考虑下列对偶问题: 求 $\Phi \in H_0^1(\Omega) \cap H^2(\Omega)$ , 使得

$$
- \nabla \cdot (\kappa \nabla \Phi) = \varepsilon_ {0}, \text {在} \Omega \text {上}. \tag {9.74}
$$

假设上述对偶问题有 $H^2$ 正则性，即存在常数 $C$ 使得

$$
\| \Phi \| _ {2} \leqslant C \| \varepsilon_ {0} \|. \tag {9.75}
$$

定理9.5 假设 $u_{h}\in V_{h}$ 是弱有限元数值解，精确解 $u\in H^{2}(\varOmega)$ 那么存在正常数 $C$ 使得

$$
\| u - u _ {0} \| \leqslant C h ^ {2} \| u \| _ {2}. \tag {9.76}
$$

证明 在方程 (9.74) 两端与 $\varepsilon_0$ 作内积, 并结合 $\sum_{K \in \mathcal{T}_h} \langle \kappa \nabla \Phi \cdot \pmb{n}, \varepsilon_b \rangle_{\partial K} = 0$ 可以得到

$$
\begin{array}{l} \left\| \varepsilon_ {0} \right\| ^ {2} = - (\nabla \cdot (\kappa \nabla \Phi), \varepsilon_ {0}) \\ = \sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla \Phi , \nabla \varepsilon_ {0}) _ {K} - \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa \nabla \Phi \cdot \boldsymbol {n}, \varepsilon_ {0} - \varepsilon_ {b} \right\rangle_ {\partial K}. \tag {9.77} \\ \end{array}
$$

在 (9.66) 中令 $u = \Phi, v = \varepsilon_h$ 得到

$$
\sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla \Phi , \nabla \varepsilon_ {0}) _ {K} = \sum_ {K \in \mathcal {T} _ {h}} (\kappa \nabla_ {w} \Phi , \nabla_ {w} \varepsilon_ {h}) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa \mathbb {Q} _ {h} \nabla \Phi \cdot \boldsymbol {n}, \varepsilon_ {0} - \varepsilon_ {b} \right\rangle_ {\partial K}. \tag {9.78}
$$

将(9.78)代入(9.77)可得

$$
\begin{array}{l} \left\| \varepsilon_ {0} \right\| ^ {2} = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} \varepsilon_ {h}, \nabla_ {w} \Phi\right) _ {K} - \sum_ {K \in \mathcal {T} _ {h}} \left\langle \kappa (\nabla \Phi - \mathbb {Q} _ {h} \nabla \Phi) \cdot \boldsymbol {n}, \varepsilon_ {0} - \varepsilon_ {b} \right\rangle_ {\partial K} \\ = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} e _ {h}, \nabla_ {w} \Phi\right) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} \left(Q _ {h} u - u\right), \nabla_ {w} \Phi\right) _ {K} - \ell_ {\Phi} (\varepsilon_ {h}) \\ = \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} e _ {h}, \nabla_ {w} Q _ {h} \Phi\right) _ {K} + \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} e _ {h}, \nabla_ {w} \left(\Phi - Q _ {h} \Phi\right)\right) _ {K} + \\ \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} \left(Q _ {h} u - u\right), \nabla_ {w} \Phi\right) _ {K} - \ell_ {\Phi} \left(\varepsilon_ {h}\right) \\ = \ell_ {u} \left(Q _ {h} \Phi\right) + \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} e _ {h}, \nabla_ {w} \left(\Phi - Q _ {h} \Phi\right)\right) _ {K} + \\ \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} \left(Q _ {h} u - u\right), \nabla_ {w} \Phi\right) _ {K} - \ell_ {\Phi} \left(\varepsilon_ {h}\right) \\ = I _ {1} + I _ {2} + I _ {3} + I _ {4}. \tag {9.79} \\ \end{array}
$$

接下来, 我们估计 (9.79) 右端的每一项. 根据 Cauchy-Schwarz 不等式, 迹不等式 (6.77), 投影算子 $Q_{h}$ 和 $\mathbb{Q}_{h}$ 的定义以及投影不等式 (9.19) 可得

$$
\begin{array}{l} I _ {1} = \left| \ell_ {u} \left(Q _ {h} \Phi\right) \right| \\ \leqslant \left| \sum_ {K \in \mathcal {T} _ {h}} \langle \boldsymbol {\kappa} (\nabla u - \mathbb {Q} _ {h} \nabla u) \cdot \boldsymbol {n}, Q _ {0} \Phi - Q _ {b} \Phi \rangle_ {\partial K} \right. \\ \leqslant C \left(\sum_ {K \in \mathcal {T} _ {h}} \| \nabla u - \mathbb {Q} _ {h} \nabla u \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {K \in \mathcal {T} _ {h}} \| Q _ {0} \Phi - Q _ {b} \Phi \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \\ \leqslant C \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} \| \nabla u - \mathbb {Q} _ {h} \nabla u \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \left(\sum_ {K \in \mathcal {T} _ {h}} h _ {K} ^ {- 1} \| Q _ {0} \Phi - \Phi \| _ {\partial K} ^ {2}\right) ^ {\frac {1}{2}} \\ \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \\ \end{array}
$$

由估计(9.69)和(9.70)可得

$$
\begin{array}{l} I _ {2} = \left| \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} e _ {h}, \nabla_ {w} \left(\Phi - Q _ {h} \Phi\right)\right) _ {K} \right| \\ \leqslant C \| e _ {h} \| \| \Phi - Q _ {h} \Phi \| \\ \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \\ \end{array}
$$

为了估计 $I_{3}$ ，我们定义 $R_{h}$ 为从 $[L^2 (K)]^d$ 到 $[P_1(K)]^d$ 的 $L^2$ 投影算子.根据弱梯度的定义(9.51）可得

$$
\begin{array}{l} \left(\kappa \nabla_ {w} \left(Q _ {h} u - u\right), R _ {h} \nabla \Phi\right) _ {K} = - \left(Q _ {0} u - u, \nabla \cdot (\kappa R _ {h} \nabla \Phi)\right) _ {K} + \left\langle Q _ {b} u - u, \kappa R _ {h} \nabla \Phi \cdot n \right\rangle_ {\partial K} \\ = 0. \\ \end{array}
$$

结合上式和估计(9.69)，投影算子 $R_{h}$ 的定义以及投影不等式(9.19)，我们有

$$
\begin{array}{l} I _ {3} = \left| \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} \left(Q _ {h} u - u\right), \nabla_ {w} \Phi\right) _ {K} \right| \\ = \left| \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} \left(Q _ {h} u - u\right), \nabla_ {w} \Phi - R _ {h} \nabla \Phi\right) _ {K} \right| \\ = \left| \sum_ {K \in \mathcal {T} _ {h}} \left(\kappa \nabla_ {w} \left(Q _ {h} u - u\right), \nabla \Phi - R _ {h} \nabla \Phi\right) _ {K} \right| \\ \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \\ \end{array}
$$

根据估计 (9.68)一(9.70)可得

$$
\begin{array}{l} I _ {4} = \left| \ell_ {\Phi} \left(\varepsilon_ {h}\right) \right| \\ \leqslant C h \| \Phi \| _ {2} \| \varepsilon_ {h} \| \\ \leqslant C h \| \Phi \| _ {2} (\| e _ {h} \| + \| u - Q _ {h} u \|) \\ \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}. \\ \end{array}
$$

结合上述所有估计与 (9.79) 可得

$$
\left\| \varepsilon_ {0} \right\| ^ {2} \leqslant C h ^ {2} \| u \| _ {2} \| \Phi \| _ {2}.
$$

再根据 $H^2$ 正则性估计（9.75）可得

$$
\| \varepsilon_ {0} \| \leqslant C h ^ {2} \| u \| _ {2}.
$$

最后由三角不等式和投影不等式 (9.18), 我们有

$$
\left\| e _ {0} \right\| \leqslant \left\| \varepsilon_ {0} \right\| + \left\| u - Q _ {0} u \right\| \leqslant C h ^ {2} \| u \| _ {2}.
$$

定理得证

# 第10章

# 有限元多重网格法

多重网格法是求解离散椭圆型问题的一种最优复杂性的算法 [17, 18, 22]. 可以用 $O(N)$ 的计算量得到与有限元解相同精度的近似解, 这里 $N$ 是有限元方程组的未知数个数.

多重网格算法的思想可以用两句话描述：在当前网格磨光；在粗网格校正。磨光步可以减少误差的高频部分；校正步利用粗网格求得的误差修正近似解，提高精度。

# 10.1 模型问题

设 $\Omega \subset \mathbb{R}^d$ ( $d = 2,3$ ) 是凸多边形或凸多面体区域. 设

$$
a (u, v) = \int_ {\Omega} \kappa \nabla u \cdot \nabla v \mathrm {d} x, \tag {10.1}
$$

其中系数 $\kappa$ 满足： $\kappa \in C^{1}(\bar{\Omega})^{d\times d}$ 且存在常数 $a_0 > 0$ 使得矩阵 $\kappa$ 满足

$$
\left(\kappa (\boldsymbol {x}) \boldsymbol {\xi}, \boldsymbol {\xi}\right) \geqslant a _ {0} (\boldsymbol {\xi}, \boldsymbol {\xi}), \quad \forall \boldsymbol {\xi} \in \mathbb {R} ^ {d}, \boldsymbol {x} \in \Omega .
$$

给定 $f\in L^{2}(\varOmega)$ ，考虑Dirichlet边值问题：求 $u\in V = H_0^1 (\varOmega)$ 使得

$$
a (u, v) = (f, v), \quad \forall v \in V. \tag {10.2}
$$

设 $\mathcal{M}_k, k = 1,2,\dots$ 是 $\Omega$ 的一列嵌套的三角剖分， $\mathcal{M}_k$ 由 $\mathcal{M}_{k-1}$ ( $k > 1$ ) 通过一致加密得到. 设 $V_k \subset H_0^1(\Omega)$ 是 $\mathcal{M}_k$ 上的连续线性有限元空间. 显然 $V_{k-1} \subset V_k$ ( $k > 1$ ). 空间 $V_k$ 上的有限元离散为：求 $u_k \in V_k$ 使得

$$
a \left(u _ {k}, v _ {k}\right) = (f, v _ {k}), \quad \forall v _ {k} \in V _ {k}. \tag {10.3}
$$

我们引入到 $V_{k}$ 的 $L^2$ 和 $H^{1}$ 投影算子

$$
\left(Q _ {k} \varphi , v _ {k}\right) = (\varphi , v _ {k}), \quad a \left(P _ {k} \psi , v _ {k}\right) = a \left(\psi , v _ {k}\right), \quad \forall v _ {k} \in V _ {k},
$$

其中 $\varphi \in L^{2}(\Omega),\psi \in H_{0}^{1}(\Omega)$ 显然 $u_{k} = P_{k}u.$ 记 $h_k = \max_{K\in \mathcal{M}_k}h_K,\| \cdot \| _A = a(\cdot ,\cdot)^{\frac{1}{2}}.$ 应用Aubin-Nitsche 技巧可得

$$
\left\| \left(I - P _ {k}\right) v \right\| _ {L ^ {2} (\Omega)} \leqslant C h _ {k} \| \left(I - P _ {k}\right) v \| _ {A}, \quad \forall v \in H _ {0} ^ {1} (\Omega). \tag {10.4}
$$

# 10.2 经典迭代法

# 10.2.1 矩阵形式和算子形式

记 $\{\phi_k^i,i = 1,2,\dots ,n_k\}$ 为 $V_{k}$ 的节点基.任给 $v_{k} = \sum_{i = 1}^{n_{k}}v_{k}^{i}\phi_{k}^{i}\in V_{k}$ 和 $g_{k}\in V_{k}$ ，如下定义 $\widetilde{v}_k = (\widetilde{v}_k^i)_{n_k\times 1},\widetilde{\widetilde{g}}_k = (\widetilde{\widetilde{g}}_k^i)_{n_k\times 1}\in \mathbb{R}^{n_k}$

$$
\widetilde {v} _ {k} ^ {i} = v _ {k} ^ {i}, \quad \widetilde {\widetilde {g}} _ {k} ^ {i} = \left(g _ {k}, \phi_ {k} ^ {i}\right), \quad i = 1, 2, \dots , n _ {k}. \tag {10.5}
$$

记 $\widetilde{A}_k = \left[a\left(\phi_k^j,\phi_k^i\right)\right]_{i,j = 1}^{n_k}$ 为刚度矩阵, 则 (10.3) 的矩阵形式为

$$
\widetilde {A} _ {k} \widetilde {u} _ {k} = \widetilde {\widetilde {f}} _ {k}, \quad \text {其 中} f _ {k} := Q _ {k} f. \tag {10.6}
$$

求解 (10.6) 的线性迭代法的一般形式为: 给定初值 $\widetilde{u}^{(0)} \in \mathbb{R}^{n_k}$ ,

$$
\widetilde {u} ^ {(n + 1)} = \widetilde {u} ^ {(n)} + \widetilde {R} _ {k} \left(\widetilde {\widetilde {f}} _ {k} - \widetilde {A} _ {k} \widetilde {u} ^ {(n)}\right), \quad n = 0, 1, 2, \dots . \tag {10.7}
$$

矩阵 $\widetilde{R}_k$ 称为迭代子.我们知道，(10.7)收敛的充要条件是谱半径 $\rho \left(I - \widetilde{R}_k\widetilde{A}_k\right) < 1.$

注意到 $\widetilde{A}_k$ 是对称正定的, 把它分解 $\widetilde{A}_k = \widetilde{D} + \widetilde{L} + \widetilde{L}^{\mathrm{T}}$ , 其中 $\widetilde{D}$ 和 $\widetilde{L}$ 分别为 $\widetilde{A}_k$ 的对角线所组成的对角矩阵和 $\widetilde{A}_k$ 的下三角形矩阵. 我们回忆以下几种迭代法:

$$
\widetilde {R} _ {k} = \left\{ \begin{array}{l l} \omega I, & \text {R i c h a r d s o n}, \\ \omega \widetilde {D} ^ {- 1}, & \text {阻 尼 J a c o b i}, \\ \left(\widetilde {D} + \widetilde {L}\right) ^ {- 1}, & \text {G a u s s - S e i d e l}. \end{array} \right. \tag {10.8}
$$

引理10.1 (i)Richardson迭代法收敛的充要条件是 $0 < \omega < \frac{2}{\rho(\widetilde{A}_k)}$

(ii) 阻尼 Jacobi 迭代法收敛的充要条件是 $0 < \omega < \frac{2}{\rho(\widetilde{D}^{-1}\widetilde{A}_k)}$   
(iii) Gauss-Seidel 迭代法恒收敛.

定义算子 $A_{k}:V_{k}\to V_{k}$

$$
\left(A _ {k} w _ {k}, v _ {k}\right) = a \left(w _ {k}, v _ {k}\right), \quad \forall v _ {k} \in V _ {k}.
$$

有限元格式 (10.3) 的算子形式为

$$
A _ {k} u _ {k} = f _ {k}. \tag {10.9}
$$

易知

$$
\widetilde {\widetilde {A _ {k} v _ {k}}} = \widetilde {A} _ {k} \widetilde {v} _ {k}, \quad \forall v _ {k} \in V _ {k}. \tag {10.10}
$$

从而，(10.9)两边同时取 $\tilde{\mathbf{\Gamma}}$ 运算即可得到(10.6).如果定义线性算子 $R_{k}:V_{k}\mapsto V_{k}$ 为

$$
R _ {k} g = \sum_ {i, j = 1} ^ {n _ {k}} \left(\widetilde {R} _ {k}\right) _ {i j} \left(g, \phi_ {k} ^ {j}\right) \phi_ {k} ^ {i}, \tag {10.11}
$$

则 $\widetilde{R_kg} = \widetilde{R_k}\widetilde{\widetilde{g}}$ ，从而求解矩阵方程（10.6）的迭代法（10.7）等价于下面求解算子方程(10.9）的迭代法：给定初值 $u^{(0)}\in V_k$

$$
u ^ {(n + 1)} = u ^ {(n)} + R _ {k} \left(f _ {k} - A _ {k} u ^ {(n)}\right), \quad n = 0, 1, 2, \dots , \tag {10.12}
$$

也就是，(10.12)两边同时取~运算即可得到(10.7).误差传播算子为 $I - R_{k}A_{k}$

引理10.2 设 $P_{k}^{i}$ 为到 $\{\phi_k^i\}$ 张成的空间的投影：

$$
a \left(P _ {k} ^ {i} w _ {k}, \phi_ {k} ^ {i}\right) = a \left(w _ {k}, \phi_ {k} ^ {i}\right), \quad \forall w _ {k} \in V _ {k}. \tag {10.13}
$$

则几种迭代法的等价算子形式（见(10.12)）的迭代子满足

$$
R _ {k} g = \left\{ \begin{array}{l l} \omega \sum_ {i = 1} ^ {n _ {k}} \left(g, \phi_ {k} ^ {i}\right) \phi_ {k} ^ {i}, & \text {R i c h a r d s o n ,} \\ \omega \sum_ {i = 1} ^ {n _ {k}} P _ {k} ^ {i} A _ {k} ^ {- 1} g, & \text {阻 尼 J a c o b i ,} \\ (I - E _ {k}) A _ {k} ^ {- 1} g, & \text {G a u s s - S e i d e l .} \end{array} \right. \tag {10.14}
$$

其中 $E_{k} = (I - P_{k}^{n_{k}})\dots (I - P_{k}^{2})(I - P_{k}^{1})$

证明 以阻尼 Jacobi 方法为例, 由 (10.13) 知

$$
P _ {k} ^ {i} w _ {k} = \frac {a (w _ {k} , \phi_ {k} ^ {i})}{a (\phi_ {k} ^ {i} , \phi_ {k} ^ {i})} \phi_ {k} ^ {i}, \quad i = 1, 2, \dots , n _ {k}.
$$

阻尼Jacobi法的迭代子为

$$
\widetilde {R} _ {k} = \omega \widetilde {D} ^ {- 1} = \omega \operatorname {d i a g} \left(a \left(\phi_ {k} ^ {1}, \phi_ {k} ^ {1}\right) ^ {- 1}, a \left(\phi_ {k} ^ {2}, \phi_ {k} ^ {2}\right) ^ {- 1}, \dots , a \left(\phi_ {k} ^ {n _ {k}}, \phi_ {k} ^ {n _ {k}}\right) ^ {- 1}\right).
$$

由 (10.11),

$$
R _ {k} g = \omega \sum_ {i = 1} ^ {n _ {k}} \frac {\left(g , \phi_ {k} ^ {i}\right)}{a \left(\phi_ {k} ^ {i} , \phi_ {k} ^ {i}\right)} \phi_ {k} ^ {i} = \omega \sum_ {i = 1} ^ {n _ {k}} \frac {a \left(A _ {k} ^ {- 1} g , \phi_ {k} ^ {i}\right)}{a \left(\phi_ {k} ^ {i} , \phi_ {k} ^ {i}\right)} \phi_ {k} ^ {i} = \omega \sum_ {i = 1} ^ {n _ {k}} P _ {k} ^ {i} A _ {k} ^ {- 1} g, \quad \forall g \in V _ {k},
$$

得证.

# 10.2.2 磨光性质

众所周知, 当 $n_k$ 大时 (即网格密时), (10.8) 中的几个迭代法不能十分有效地求解有限元方程组 (10.6). 为了使得误差下降一半, 需要 $O\left(h_k^{-2}\right)$ 步迭代, 即需要 $O\left(n_k^{1 + \frac{2}{d}}\right)$ 的计算量. 但是它们都有一个重要的所谓“磨光性质”. 例如考虑 (10.6) 的Richardson迭代 $\omega = \frac{1}{\rho(\widetilde{A}_k)}$ :

$$
\widetilde {u} ^ {(n + 1)} = \widetilde {u} ^ {(n)} + \frac {1}{\rho (\widetilde {A} _ {k})} \left(\widetilde {\widetilde {f}} _ {k} - \widetilde {A} _ {k} \widetilde {u} ^ {(n)}\right), \quad n = 0, 1, 2, \dots .
$$

由于 $\widetilde{A}_k$ 是对称正定的，它存在 $n_k$ 个实特征值 $0 < \mu_1 \leqslant \mu_2 \leqslant \dots \leqslant \mu_{n_k}$ ，且可取对应的特征向量为 $\widetilde{\phi}_1, \widetilde{\phi}_2, \dots, \widetilde{\phi}_{n_k}$ 满足 $\left(\widetilde{\phi}_i, \widetilde{\phi}_j\right) = \delta_{ij}$ . 既然刚度矩阵对应椭圆型微分算子的离散，一般说来，特征值越大，其对应的特征向量振荡得越厉害. 设 $\widetilde{u}_k - \widetilde{u}^0 = \sum_{i=1}^{n_k} \alpha_i \widetilde{\phi}_i$ ，则

$$
\widetilde {u} _ {k} - \widetilde {u} ^ {(n)} = \sum_ {i} \alpha_ {i} \left(1 - \frac {\mu_ {i}}{\mu_ {n _ {k}}}\right) ^ {n} \widetilde {\phi} _ {i}.
$$

显然，如果 $\mu_{i}$ 接近 $\mu_{n_k}$ ，那么当 $n\to \infty$ 时， $\left(1 - \frac{\mu_i}{\mu_{n_k}}\right)^n$ 以很快的速度趋于0.这意味着误差的高频分量衰减很快

下面给出一个数值例子演示一下 Gauss-Seidel 方法的磨光性质. 考虑单位正方形上的 Poisson 问题 $-\Delta u = 1, x \in \Omega, u = 0, x \in \partial \Omega.$ 取一致三角剖分. 图 10.1 显示误差的高频部分在 Gauss-Seidel 迭代中衰减很快. Gauss-Seidel 迭代收敛慢是因为对误差的低频部分衰减得不好.

![](images/df2d7d1827a58d2d6af311351f68e0660327853a8937547530c52b8b3f8e3acb.jpg)

![](images/9dd5165872a7ced565bb271fb84dde2c48091d52e86ab51cfec0434186dabf0c.jpg)

![](images/4b1d5b580a02098bc445ccb6cfd16768048e1895d81544b34db363aedda49912.jpg)

![](images/6ecabcaa1371982b40f7d9cdd950cde6cc2b57f7e92be86cea3f95ffb7190c8e.jpg)  
图10.1 Gauss-Seidel迭代0,3,9和200次以后的误差，未知数个数为2113

对于上面的模型问题, Brandt 利用“局部模态分析”的办法得知阻尼 Jacobi 方法当 $\omega = \frac{4}{5}$ 具有最好的磨光性质; Gauss-Seidel 方法的磨光效果比阻尼 Jacobi 方法好; 红黑

序的 Gauss-Seidel 方法比字典序的要好. 我们知道 Jacobi 方法和红黑序的 Gauss-Seidel 适合于并行计算.

下面引理刻画了阻尼Jacobi方法的磨光性质，由于篇幅所限，略去证明.

引理10.3 假设 $\omega > 0$ 足够小. 记 $K_{k} = I - R_{k}A_{k}$ , 则存在常数 $\alpha > 0$ , 使得对自然数 $m > 0$ ,

$$
\left\| \left(\boldsymbol {I} - P _ {k - 1}\right) K _ {k} ^ {m} v \right\| _ {A} ^ {2} \leqslant \frac {\alpha}{m} \left(\| v \| _ {A} ^ {2} - \| K _ {k} ^ {m} v \| _ {A} ^ {2}\right), \quad \forall v \in V _ {k}.
$$

此引理说明，阻尼 Jacobi 迭代后，近似解的误差可以被上一层的较粗网格上的有限元函数很好地逼近。正因为经典迭代法的磨光效应，我们把其迭代子也称为光滑子。

# 10.3 多重网格 $V$ 循环算法

多重网格法的基本思想是在当前网格磨光, 在粗网格校正. 设 $R_{k}:V_{k}\to V_{k}$ 为一光滑子, $R_{k}^{t}$ 是 $R_{k}$ 关于 $L^2$ 内积 $(\cdot ,\cdot)$ 的伴随算子. 求解 (10.9) 的多重网格 $V$ 循环算法可写为一般线性迭代法的形式: 给定 $u^{(0)}\in V_k$

$$
u ^ {(n + 1)} = u ^ {(n)} + \mathbb {B} _ {k} \left(f _ {k} - A _ {k} u ^ {(n)}\right), \quad n = 0, 1, 2, \dots , \tag {10.15}
$$

关键是如何定义迭代子 $\mathbb{B}_k$ 。显然，只需对任意的 $g \in V_k$ 定义 $\mathbb{B}_kg$ 。类比（10.15）（若 $u^{(0)} = 0$ ，则 $u^{(1)} = \mathbb{B}_kf_k$ ），相当于求解方程 $A_ky = g$ ，以零为初值按多重网格方法迭代一次得到的近似解即为 $\mathbb{B}_kg$ 。

为此我们先用经典迭代法从 $y^0 = 0$ 出发迭代几次做磨光，比如迭代 $m$ 次得到 $y^{m}$

记误差为 $e = y - y^{m}$ 显然 $\mathcal{E}$ 满足 $A_{k}e = r\coloneqq g - A_{k}y^{m}$ ，正好是右端为残量 $\pmb{r}$ 的有限元解，等价地写为：求 $e\in V_k$ 使得

$$
a (e, v) = (r, v), \quad \forall v \in V _ {k}. \tag {10.16}
$$

如果在细网格空间 $V_{k}$ 上把误差 $e$ 精确求出来，用它来校正近似解 $y^{m}$ ，即计算 $y^{m} + e$ 正好就能得到精确解 $y$ 。但这样做代价太大，不亚于直接求解原方程 $A_{k}y = g$ 。既然磨光后误差 $e$ 变得光滑了，用粗网格函数就可以逼近得比较好了，所以自然地，我们可以把(10.16)限制到粗网格空间 $V_{k-1}$ 来近似求解误差，即：求 $e \in V_{k-1}$ 使得

$$
a (e, v) = (r, v) = \left(Q _ {k - 1} r, v\right), \quad \forall v \in V _ {k - 1}.
$$

其等价的算子形式为

$$
A _ {k - 1} e = Q _ {k - 1} r. \tag {10.17}
$$

求解上面的粗网格方程得到 $e \in V_{k-1}$ 应该是上面细网格误差的很好近似，而且显著减少了计算量。用这个 $e$ 校正 $y^{m}$ 得到 $y^{m} + e$ 应该也是 $y$ 的一个好的近似。如果精确求解粗网格误差方程 (10.17)，就得到“两重网格方法”。当然，可以继续用多重网格方法求解(10.17)，即以零为初值（误差 $e$ 一般很小），迭代一次得到 $\mathbb{B}_{k-1}Q_{k-1}r$ 作为误差 $e$ 的近似。这样递归下去直到最粗网格精确求解，就得到了多重网格方法。

总结上述过程就可以得到多重网格迭代子的递归算法10.1，包含三个步骤：细网格前磨光，粗网格校正，后磨光。前两个步骤如前所述，后磨光主要是为了对称性。

# 算法10.1 $V$ 循环迭代子

对 $k = 1$ ，定义 $\mathbb{B}_1 = A_1^{-1}$ .假设 $\mathbb{B}_{k - 1}:V_{k - 1}\mapsto V_{k - 1}$ 已定义.迭代子 $\mathbb{B}_k:V_k\mapsto V_k$ 定义如下.任给 $g\in V_k$

1. 前磨光: 取 $y^0 = 0 \in V_k$ , 对 $j = 1,2,\dots,m$

$$
y ^ {j} = y ^ {j - 1} + R _ {k} \left(g - A _ {k} y ^ {j - 1}\right).
$$

2. 粗网格校正: $e = \mathbb{B}_{k-1}Q_{k-1}\left(g - A_ky^m\right), y^{m+1} = y^m + e,$

3. 后磨光: 对 $j = m + 2, m + 3, \dots, 2m + 1$ ,

$$
y ^ {j} = y ^ {j - 1} + R _ {k} ^ {t} \left(g - A _ {k} y ^ {j - 1}\right).
$$

定义 $\mathbb{B}_k g = y^{2m + 1}$

顺便说明一下，多重网格迭代子也可以被用来作为预条件子处理一下有限元方程改善其条件数，再用其他迭代法求解，比如共轭斜量法，得到多重网格预处理的共轭斜量(MG-PCG)法.

为了证明收敛性，我们先给出多重网格算法误差传播算子的递推关系.记 $y = A_k^{-1}g$ 则

$$
y - y ^ {2 m + 1} = \left(\boldsymbol {I} - R _ {k} ^ {t} A _ {k}\right) ^ {m} \left(\boldsymbol {I} - \mathbb {B} _ {k - 1} Q _ {k - 1} A _ {k}\right) \left(\boldsymbol {I} - R _ {k} A _ {k}\right) ^ {m} \left(y - y ^ {0}\right).
$$

记 $K_{k} = I - R_{k}A_{k}$ ， $K_{k}^{*}$ 为 $K_{k}$ 关于内积 $a(\cdot ,\cdot)$ 的伴随算子，则 $K_{k}^{*} = I - R_{k}^{t}A_{k}$ 且

$$
\boldsymbol {I} - \mathbb {B} _ {k} A _ {k} = K _ {k} ^ {* m} \left(\boldsymbol {I} - \mathbb {B} _ {k - 1} Q _ {k - 1} A _ {k}\right) K _ {k} ^ {m}.
$$

另外，对任何 $v_{k}\in V_{k},w_{k - 1}\in V_{k - 1}$ ，有

$$
\begin{array}{l} \left(Q _ {k - 1} A _ {k} v _ {k}, w _ {k - 1}\right) = \left(A _ {k} v _ {k}, w _ {k - 1}\right) = a \left(v _ {k}, w _ {k - 1}\right) \\ = a \left(P _ {k - 1} v _ {k}, w _ {k - 1}\right) = \left(A _ {k - 1} P _ {k - 1} v _ {k}, w _ {k - 1}\right), \\ \end{array}
$$

从而 $Q_{k - 1}A_{k} = A_{k - 1}P_{k - 1}$ 因此我们得到下面的递推关系：

引理10.4 在 $V_{k}$ 中成立

$$
\boldsymbol {I} - \mathbb {B} _ {k} A _ {k} = K _ {k} ^ {* m} \left(\left(\boldsymbol {I} - P _ {k - 1}\right) + \left(\boldsymbol {I} - \mathbb {B} _ {k - 1} A _ {k - 1}\right) P _ {k - 1}\right) K _ {k} ^ {m}.
$$

下面仅对采用阻尼Jacobi磨光的多重网格 $V$ 循环算法, 给出收敛性分析. 此时

$$
R _ {k} = \omega \sum_ {i = 1} ^ {n _ {k}} P _ {k} ^ {i} A _ {k} ^ {- 1}.
$$

下面引理的证明留作习题

引理10.5 对采用阻尼Jacobi磨光的多重网格 $V$ 循环算法有

$$
R _ {k} ^ {t} = R _ {k}, \quad K _ {k} ^ {*} = K _ {k}, \quad \left(R _ {k} v _ {k}, v _ {k}\right) \geqslant 0, \quad \forall v _ {k} \in V _ {k};
$$

$$
B _ {k} ^ {t} = B _ {k}, \quad (I - \mathbb {B} _ {k} A _ {k}) ^ {*} = (I - \mathbb {B} _ {k} A _ {k}).
$$

定理10.1 在引理10.3的条件下，有如下的关于带阻尼Jacobi磨光的多重网格 $V$ 循环算法的收敛性估计：

$$
\| I - \mathbb {B} _ {k} A _ {k} \| _ {A} \leqslant \delta := \frac {\alpha}{\alpha + m}, \tag {10.18}
$$

其中 $\alpha > 0$ 与网格 $\mathcal{M}_k$ 及磨光次数 $m \geqslant 1$ 无关，

$$
\| I - \mathbb {B} _ {k} A _ {k} \| _ {A} = \sup  _ {0 \neq v \in V _ {k}} \frac {\| (I - \mathbb {B} _ {k} A _ {k}) v \| _ {A}}{\| v \| _ {A}} = \sup  _ {0 \neq v \in V _ {k}} \frac {a ((I - \mathbb {B} _ {k} A _ {k}) v , v)}{\| v \| _ {A} ^ {2}}.
$$

证明 我们用数学归纳法证明：

$$
0 \leqslant a \left(\left(I - \mathbb {B} _ {k} A _ {k}\right) v, v\right) \leqslant \delta a (v, v), \quad \forall v \in V _ {k}. \tag {10.19}
$$

既然 $\mathbb{B}_1 = A_1^{-1}$ , $k = 1$ 时显然成立. 假设 $k - 1$ 时成立, 下面证明 $k$ 时成立. 由引理 10.4-10.5, 对任意 $v \in V_k$ ,

$$
\begin{array}{l} a \left(\left(I - \mathbb {B} _ {k} A _ {k}\right) v, v\right) = a \left(K _ {k} ^ {m} \left(I - P _ {k - 1}\right) K _ {k} ^ {m} v, v\right) + \\ a \left(K _ {k} ^ {m} \left(I - \mathbb {B} _ {k - 1} A _ {k - 1}\right) P _ {k - 1} K _ {k} ^ {m} v, v\right) \\ = a \left(\left(I - P _ {k - 1}\right) K _ {k} ^ {m} v, \left(I - P _ {k - 1}\right) K _ {k} ^ {m} v\right) + \\ a \left(\left(I - \mathbb {B} _ {k - 1} A _ {k - 1}\right) P _ {k - 1} K _ {k} ^ {m} v, P _ {k - 1} K _ {k} ^ {m} v\right) \\ \geqslant 0. \\ \end{array}
$$

另外，由引理10.3

$$
\begin{array}{l} a \left(\left(I - \mathbb {B} _ {k} A _ {k}\right) v, v\right) \leqslant \| \left(I - P _ {k - 1}\right) K _ {k} ^ {m} v \| _ {A} ^ {2} + \delta \| P _ {k - 1} K _ {k} ^ {m} v \| _ {A} ^ {2} \\ = (1 - \delta) \| (I - P _ {k - 1}) K _ {k} ^ {m} v \| _ {A} ^ {2} + \delta \| K _ {k} ^ {m} v \| _ {A} ^ {2} \\ \leqslant (1 - \delta) \frac {\alpha}{m} \left(\| v \| _ {A} ^ {2} - \| K _ {k} ^ {m} v \| _ {A} ^ {2}\right) + \delta \| K _ {k} ^ {m} v \| _ {A} ^ {2} \\ = \delta \| v \| _ {A} ^ {2}. \\ \end{array}
$$

证毕.

我们指出, 定理 10.1 的结论对带 Gauss-Sediel 磨光的多重网格 $V$ 循环算法也成立.

例10.1 考虑单位正方形 $\Omega$ 上的 Poisson 问题 $-\Delta u = 1, x \in \Omega, u|_{\partial \Omega} = 0$ 的线性有限元离散，采用三角剖分。初始剖分有 4 个三角形组成。我们用 $V$ 循环算法 (10.15) 求解有限元方程组。初值为零，采用 Gauss-Seidel 磨光，磨光次数 $m = 2$ ，终止条件为

$$
\left\| \widetilde {\widetilde {f}} _ {k} - \widetilde {A} _ {k} \widetilde {u} _ {k} ^ {(n)} \right\| _ {\infty} / \left\| \widetilde {\widetilde {f}} _ {k} - \widetilde {A} _ {k} \widetilde {u} _ {k} ^ {(0)} \right\| _ {\infty} <   1 0 ^ {- 6}.
$$

表10.1给出了1—10次一致加密的多重网格法的迭代次数．最后一层网格含4194304个三角形，有2095105个内部节点.

表 10.1 例 10.1: 未知数个数 $N$ 和多重网格法的迭代次数 $l$   

<table><tr><td>N</td><td>5</td><td>25</td><td>113</td><td>481</td><td>1985</td><td>8065</td><td>32513</td><td>130561</td><td>523265</td><td>2095105</td></tr><tr><td>l</td><td>3</td><td>6</td><td>6</td><td>7</td><td>7</td><td>7</td><td>7</td><td>7</td><td>7</td><td>7</td></tr></table>

# 10.4 完全多重网格法和工作量估计

我们知道多重网格 $V$ 循环算法的收敛速度是最优阶的, 即与加密次数 $k$ 无关. 在这一节, 我们将进一步证明, 多重网格 $V$ 循环算法的每步迭代的计算量也是最优阶的, 即 $O(n_{k})$ . 但是, 要达到 $O(h_{k})$ 的误差, 多重网格 $V$ 循环算法需要 $O\left(\log \frac{1}{h_k}\right) = O(\log n_k)$ 次迭代, 从而需要 $O(n_{k}\log n_{k})$ 的计算量. 这不是最优阶的. 在本节, 我们还将引入“完全多重网格方法”, 仅用 $O(n_{k})$ 的计算量, 就可以得到精度为 $O(h_{k})$ 的近似解, 是真正的最优阶方法.

由第3章的收敛性理论，我们知道第 $k$ 层网格上的有限元解 $u_{k}$ 满足如下误差估计：

$$
\left\| u - u _ {k} \right\| _ {A} \leqslant c _ {1} h _ {k}, \quad k \geqslant 1, \tag {10.20}
$$

其中 $c_{1} > 0$ 是与 $k$ 无关的常数

完全多重网格法的设计基于以下事实: $u_{k-1} \in V_{k-1} \subset V_k$ 是 $u_k \in V_k$ 的近似, 因此可以被用于求解 $u_k$ 的多重网格迭代法的初值. 具体算法见算法 10.2.

# 算法10.2 完全多重网格

给定正整数 $l$

对 $k = 1$ ，计算 $\hat{u}_1 = A_1^{-1}f_1$

对 $k\geqslant 2$ ，令 $\hat{u}_k = \hat{u}_{k - 1}$ ，作 $l$ 次多重网格迭代： $\hat{u}_k\gets \hat{u}_k + \mathbb{B}_k(f_k - A_k\hat{u}_k)$

记 $\tilde{h}_k = \max_{K\in \mathcal{M}_k}|K|^{\frac{1}{d}}$ 易知存在 $p > 1$ 使得 $\tilde{h}_k = \frac{\tilde{h}_{k - 1}}{p}$ 显然 $\tilde{h}_k$ 与 $h_k$ 等价，即存在

仅依赖于网格的正则性的正数 $c_{2}$ 和 $c_{3}$ 使得 $c_{2}\tilde{h}_{k} \leqslant h_{k} \leqslant c_{3}\tilde{h}_{k}$ .

下面定理说明，完全多重网格法得到的近似解的 $H^{1}$ 误差与有限元解具有同样的收敛阶.

定理10.2 假设（10.18）成立且 $\delta^l < \frac{1}{p}$ 则

$$
\left\| u _ {k} - \hat {u} _ {k} \right\| _ {A} \leqslant \frac {c _ {3} p \delta^ {l}}{c _ {2} \left(1 - p \delta^ {l}\right)} c _ {1} h _ {k}, \quad k \geqslant 1.
$$

证明 由(10.18)，

$$
\begin{array}{l} \left\| u _ {k} - \hat {u} _ {k} \right\| _ {A} \leqslant \delta^ {l} \left\| u _ {k} - \hat {u} _ {k - 1} \right\| _ {A} \leqslant \delta^ {l} \left(\left\| u _ {k} - u _ {k - 1} \right\| _ {A} + \left\| u _ {k - 1} - \hat {u} _ {k - 1} \right\| _ {A}\right) \\ \leqslant \delta^ {l} \left(\| u - u _ {k - 1} \| _ {A} + \| u _ {k - 1} - \hat {u} _ {k - 1} \| _ {A}\right) \\ \leqslant \delta^ {l} \| u _ {k - 1} - \hat {u} _ {k - 1} \| _ {A} + c _ {1} \delta^ {l} h _ {k - 1}. \\ \end{array}
$$

注意到 $\| u_1 - \hat{u}_1\| _A = 0$ ，可得

$$
\begin{array}{l} \left\| u _ {k} - \hat {u} _ {k} \right\| _ {A} \leqslant c _ {1} \sum_ {n = 1} ^ {k - 1} \left(\delta^ {l}\right) ^ {n} h _ {k - n} \leqslant c _ {1} c _ {3} \sum_ {n = 1} ^ {k - 1} \left(\delta^ {l}\right) ^ {n} \tilde {h} _ {k - n} \\ \leqslant c _ {1} c _ {3} \tilde {h} _ {k} \sum_ {n = 1} ^ {k - 1} \left(p \delta^ {l}\right) ^ {n} \leqslant \frac {c _ {1} c _ {3}}{c _ {2}} \frac {p \delta^ {l}}{1 - p \delta^ {l}} h _ {k}. \\ \end{array}
$$

证毕.

下面我们转到计算量估计. 显然

$$
n _ {k} = \dim V _ {k} \sim \frac {1}{h _ {k} ^ {d}} \sim \frac {1}{\tilde {h} _ {k} ^ {d}} \sim \left(p ^ {d}\right) ^ {k}. \tag {10.21}
$$

定理10.3 完全多重网格法的工作量是 $O(n_{k})$

证明 令 $W_{k}$ 表示第 $k$ 层 $V$ 循环迭代的计算量. 易知

$$
W _ {k} \leqslant C m n _ {k} + W _ {k - 1}.
$$

因此

$$
W _ {k} \leqslant C m \left(n _ {1} + n _ {2} + \dots + n _ {k}\right) \leqslant C n _ {k}.
$$

记 $\hat{W}_k$ 为完全多重网格方法中求得 $\hat{u}_k$ 的计算量, 则

$$
\hat {W} _ {k} \leqslant \hat {W} _ {k - 1} + l W _ {k} \leqslant \hat {W} _ {k - 1} + C n _ {k}.
$$

因此

$$
\hat {W} _ {k} \leqslant C (n _ {1} + n _ {2} + \dots + n _ {k}) \leqslant C n _ {k}.
$$

证毕.

# 10.5 多重网格 $V$ 循环算法的矩阵形式

为了便于程序实现，本节利用算符~和~将多重网格 $V$ 循环算法及完全多重网格算法的算子形式改写为矩阵形式.

记 $\{\phi_k^1,\phi_k^2,\dots ,\phi_k^{n_k}\}$ 为 $V_{k}$ 的节点基，我们定义所谓的“延拓矩阵” $I_{k - 1}^{k}\in \mathbb{R}^{n_{k}\times n_{k - 1}}$ 如下：

$$
\phi_ {k - 1} ^ {j} = \sum_ {i = 1} ^ {n _ {k}} \left(I _ {k - 1} ^ {k}\right) _ {i j} \phi_ {k} ^ {i}. \tag {10.22}
$$

由 $\widetilde{v}_k$ 和 $\widetilde{\widetilde{v}}_k$ 的定义 (10.5) 得

$$
\widetilde {v} _ {k} = I _ {k - 1} ^ {k} \widetilde {v} _ {k - 1}, \quad \forall v _ {k} = v _ {k - 1}, v _ {k} \in V _ {k}, v _ {k - 1} \in V _ {k - 1}, \tag {10.23}
$$

$$
\widetilde {Q _ {k - 1} r _ {k}} = \left(I _ {k - 1} ^ {k}\right) ^ {\mathrm {T}} \widetilde {\widetilde {r}} _ {k}, \quad \forall r _ {k} \in V _ {k}. \tag {10.24}
$$

上面这两个性质表明，通过延拓矩阵 $I_{k - 1}^{k}$ 可以将粗网格解向量插值得到细网格解向量；通过 $\left(I_{k - 1}^{k}\right)^{\mathrm{T}}$ 将细网格残差向量投影到粗网格残量.

注意到 $A_{k - 1}v_{k - 1} = Q_{k - 1}A_kv_{k - 1},\forall v_{k - 1}\in V_{k - 1}$ ，我们有

$$
\widetilde {A} _ {k - 1} \widetilde {v} _ {k - 1} = \widetilde {A _ {k - 1} v _ {k - 1}} = \left(I _ {k - 1} ^ {k}\right) ^ {\mathrm {T}} \widetilde {A _ {k} v _ {k - 1}} = \left(I _ {k - 1} ^ {k}\right) ^ {\mathrm {T}} \widetilde {A} _ {k} I _ {k - 1} ^ {k} \widetilde {v} _ {k - 1},
$$

即利用延拓矩阵对细网格刚度矩阵作合同变换就可以得到粗网格刚度矩阵：

$$
\widetilde {A} _ {k - 1} = \left(I _ {k - 1} ^ {k}\right) ^ {\mathrm {T}} \widetilde {A} _ {k} I _ {k - 1} ^ {k}. \tag {10.25}
$$

对多重网格 $V$ 循环算法 (10.15) 及迭代子算法 10.1 取 $\sim$ 运算, 并利用性质 (10.10) (10.23) (10.24), 即得下面的多重网格 $V$ 循环算法的矩阵形式: 给定 $\widetilde{u}^{(0)} \in \mathbb{R}^{n_k}$ ,

$$
\widetilde {u} ^ {(n + 1)} = \widetilde {u} ^ {(n)} + \widetilde {\mathbb {B}} _ {k} \left(\widetilde {\widetilde {f}} _ {k} - \widetilde {A} _ {k} \widetilde {u} ^ {(n)}\right), \quad n = 0, 1, 2, \dots \tag {10.26}
$$

及迭代子算法矩阵形式 10.3

# 算法10.3 $V$ 循环迭代子：矩阵形式

令 $\widetilde{\mathbb{B}}_1 = \widetilde{A}_1^{-1}$ . 假设 $\widetilde{\mathbb{B}}_{k - 1} \in \mathbb{R}^{n_{k - 1} \times n_{k - 1}}$ 已定义, 则 $\widetilde{\mathbb{B}}_k \in \mathbb{R}^{n_k \times n_k}$ 定义如下: $\forall \widetilde{g} \in \mathbb{R}^{n_k}$ .

1. 前磨光: 对 $\widetilde{y}^0 = 0$ 及 $j = 1,2,\dots ,m,$

$$
\widetilde {y} ^ {j} = \widetilde {y} ^ {j - 1} + \widetilde {R} _ {k} \left(\widetilde {\widetilde {g}} - \widetilde {A} _ {k} \widetilde {y} ^ {j - 1}\right).
$$

2. 粗网格校正: $\widetilde{e} = \widetilde{\mathbb{B}}_{k-1}\left(I_{k-1}^{k}\right)^{\mathrm{T}}\left(\widetilde{g} - \widetilde{A}_{k}\widetilde{y}^{m}\right), \widetilde{y}^{m+1} = \widetilde{y}^{m} + I_{k-1}^{k}\widetilde{e}$ .   
3. 后磨光: 对于 $j = m + 2, m + 3, \dots, 2m + 1$ ,

$$
\widetilde {y} ^ {j} = \widetilde {y} ^ {j - 1} + \widetilde {R} _ {k} ^ {\mathrm {T}} \left(\widetilde {\widetilde {g}} - \widetilde {A} _ {k} \widetilde {y} ^ {j - 1}\right).
$$

定义 $\tilde{\mathbb{B}}_k\tilde{\tilde{g}} = \tilde{y}^{2m + 1}$

类似可得完全多重网格算法10.2的矩阵形式算法10.4

# 算法10.4 完全多重网格算法：矩阵形式

对 $k = 1,\widetilde{u}_1 = \widetilde{A}_1^{-1}\widetilde{\widetilde{f}}_1$

对 $k \geqslant 2$ , 令 $\widetilde{u}_k = I_{k-1}^k \widetilde{u}_{k-1}$ , 作 $l$ 次迭代 $\widetilde{u}_k \gets \widetilde{u}_k + \widetilde{\mathbb{B}}_k \left( \widetilde{\widetilde{f}_k} - \widetilde{A}_k \widetilde{u}_k \right)$ .

# 10.5.1 习题

1. 证明引理 10.2 的 Gauss-Seidel 迭代的情形  
2. 证明引理10.5   
3. 利用完全多重网格法求解如下椭圆型问题的线性有限元离散：

$$
\left\{ \begin{array}{l l} {- \nabla \cdot (a \nabla u) = 1,} & {x \in \varOmega := (- 1, 1) \times (- 1, 1),} \\ {u = 0,} & {x \in \partial \varOmega ,} \end{array} \right. a = \left\{ \begin{array}{l l} {1,} & {x _ {1} x _ {2} > 0,} \\ {1 0,} & {x _ {1} x _ {2} \leqslant 0.} \end{array} \right.
$$

探究不同的 $l$ 的选取对完全多重网格法所求近似解的精度的影响 (参见定理 10.2,与直接法的解或与 $V$ 循环多重网格解比较即可).

# 第11章

# 自适应有限元法

基于可计算的后验误差估计的自适应有限元法（参见[18,35]）提供了一种系统的局部加密或粗化网格的方法，适于求解带奇性的定常或发展型偏微分问题.本章的内容取自[18],目的是以求解二阶椭圆型问题为例，来描述自适应有限元法的基本思想

# 11.1 一个带奇性的例子

我们知道, 如果椭圆型问题的解 $Lu = f$ 具有正则性 $u \in H^{2}(\Omega)$ , 则线性有限元法可以达到最优收敛阶 $O(h)$ . 但是, 对于具有凹角的区域, 解一般不再属于 $H^{2}(\Omega)$ . 因此, 经典的拟均匀网格上的有限元法无法给出令人满意的结果, 即达不到最优阶收敛. 本章的目的是提供一种解决这一问题的方法.

首先我们给出一个例子来展示凹角带来奇性行为. 给定角度 $0 < \omega < 2\pi$ ，我们考虑扇形区域 $S_{\omega} = \{(r, \theta) : 0 < r < \infty, 0 \leqslant \theta \leqslant \omega\}$ （如图11.1）上的调和函数，即 $-\Delta u = 0$ 满足边界条件： $u|_{\Gamma_1 \cup \Gamma_2} = 0$ ，其中

$$
\Gamma_ {1} = \{(r, \theta): r > 0, \theta = 0 \},
$$

$$
\Gamma_ {2} = \{(r, \theta): r > 0, \theta = \omega \}.
$$

![](images/6c45d87a17d0a1dc8529144eab1305eb581c78578e547d02ae965266424efb96.jpg)  
图11.1 示意图

我们利用分离变量法，考虑形如 $u = r^{\alpha}\mu (\theta)$ 的解.由于在极坐标中

$$
\Delta u = \frac {\partial^ {2} u}{\partial r ^ {2}} + \frac {1}{r} \frac {\partial u}{\partial r} + \frac {1}{r ^ {2}} \frac {\partial u}{\partial \theta^ {2}},
$$

我们有

$$
\Delta u = \alpha (\alpha - 1) r ^ {\alpha - 2} \mu (\theta) + \alpha r ^ {\alpha - 2} \mu (\theta) + r ^ {\alpha - 2} \mu^ {\prime \prime} (\theta) = 0,
$$

这意味着

$$
\mu^ {\prime \prime} (\theta) + \alpha^ {2} \mu (\theta) = 0,
$$

因此 $\mu (\theta) = A\sin \alpha \theta +B\cos \alpha \theta .$ 由边界条件 $\mu (0) = \mu (\omega) = 0$ 得出 $\alpha = \frac{k\pi}{\omega}$ 和 $\mu (\theta) =$ $A\sin \left(\frac{k\pi}{\omega}\theta\right),k = 1,2,3,\dots .$ 因此边值问题： $\Delta u = 0,x\in S_{\omega};u = 0,x\in \Gamma_1\cup \Gamma_2,$ 有非平凡解

$$
u = r ^ {\alpha} \sin (\alpha \theta), \quad \alpha = \frac {\pi}{\omega}.
$$

引理11.1 如果 $\pi < \omega < 2\pi$ ，那么对于任何 $R > 0, u \notin H^{2}(S_{\omega} \cap B_{R})$

证明 注意到此时 $\alpha \in \left(\frac{1}{2}, 1\right)$ , 直接计算得

$$
\begin{array}{l} \int_ {\Omega} \left| \frac {\partial^ {2} u}{\partial r ^ {2}} \right| ^ {2} d x = \int_ {0} ^ {R} \int_ {0} ^ {\omega} | \alpha (\alpha - 1) r ^ {\alpha - 2} \sin (\alpha \theta) | ^ {2} r d r d \theta \\ = \alpha^ {2} (\alpha - 1) ^ {2} \int_ {0} ^ {\omega} | \sin^ {2} (\alpha \theta) | \mathrm {d} \theta \cdot \int_ {0} ^ {R} r ^ {2 \alpha - 3} \mathrm {d} r \\ = c \left. r ^ {2 (\alpha - 1)} \right| _ {0} ^ {R} = + \infty . \\ \end{array}
$$

得证.

例11.1 考虑如图11.2所示的 $L$ 形区域 $\Omega$ 上的Laplace方程 $\Delta u = 0$ , 取Dirichlet边界条件, 使得其精确解为 $u = r^{\frac{2}{3}}\sin \left(\frac{2}{3}\theta\right)$ .

![](images/17b0c8fb74caa6148626cef3f09259084b23af3cf827b0460ce5afa56baa1813.jpg)  
(a) $L$ 形区域和初始网格

![](images/224afb5739b3b7b57b033fe9be54f2a3d22f812de1b7f49098541b50703f7bf7.jpg)  
(b)对数-对数坐标下 $H^{1}$ 误差相对于 $\frac{1}{h}$ 的曲线点线给出参考斜率 $-\frac{2}{3}$   
图11.2

设 $\mathcal{M}_h$ 为 $\Omega$ 的一个三角剖分, $u_h$ 为 $\mathcal{M}_h$ 上的线性有限元解. 由于 $u \notin H^2(\Omega)$ , 定理6.10中关于 $u_h$ 的满阶 $H^1$ 误差估计对于这个 $L$ 形区域问题不成立. 为了展示线性有限元解 $u_h$ 的收敛速度, 我们取如图11.2(a)所示的初始网格, 并采用一致四等分加密, 图11.2(b)画出了对数-对数坐标下的 $H^1$ 误差 $\| u - u_h\|_{H^1(\Omega)}$ 关于 $\frac{1}{h}$ 的曲线. 结果表明, 对于拟均匀三角剖分上 $L$ 型区域问题的线性有限元近似, 应该由如下误差估计:

$$
\| u - u _ {h} \| _ {H ^ {1} (\Omega)} \leqslant C h ^ {\frac {2}{3}}.
$$

如果记 $N$ 为自由度数，即有限元方程组未知数的个数，则显然 $h = O(N^{-\frac{1}{2}})$ ，从而

$$
\left\| u - u _ {h} \right\| _ {H ^ {1} (\Omega)} \leqslant C N ^ {- \frac {1}{3}}. \tag {11.1}
$$

# 11.2 后验误差分析

设 $\Omega \subset R^{d}$ ( $d = 2,3$ ) 为有界的多面体区域, $\{\mathcal{M}_h\}$ 为 $\Omega$ 的一族正则三角剖分. 将网格 $\mathcal{M}_h$ 的所有内部边/面的集合记为 $\mathcal{E}_h^I$ . 设 $V_h \subset H^1(\Omega)$ 为 $\mathcal{M}_h$ 上的线性有限元空间, 记 $V_h^0 = V_h \cap H_0^1(\Omega)$ . 对任意 $K \in \mathcal{M}_h$ , 记 $h_K$ 为其直径; 对任意 $e = K_1 \cap K_2 \in \mathcal{E}_h^I$ , 记 $h_e$ 为其直径, $\Omega_e = K_1 \cup K_2$ . 需要指出的是, $\Omega$ 不一定是凸的.

考虑变分问题：求 $u\in H_0^1 (\varOmega)$ 使得

$$
(a \nabla u, \nabla v) = (f, v), \quad \forall v \in H _ {0} ^ {1} (\Omega), \tag {11.2}
$$

其中 $f\in L^{2}(\Omega)$ ，不妨假设 $a(x)$ 在 $\mathcal{M}_h$ 上是正的分片常数函数.设 $u_{h}\in V_{h}^{0}$ 为其有限元解：

$$
\left(a \nabla u _ {h}, \nabla v _ {h}\right) = \left(f, v _ {h}\right), \quad \forall v _ {h} \in V _ {h} ^ {0}. \tag {11.3}
$$

在本节中, 我们首先介绍非光滑函数的 Scott-Zhang 插值算子, 然后介绍后验误差估计, 包括上界估计和下界估计.

# 11.2.1 Scott-Zhang插值算子

我们知道，二阶椭圆型问题的弱解属于 $H^1$ 空间，但 $H^1$ 空间中的函数不一定连续(见例5.2)，所以其Lagrange插值不一定有定义.本小节介绍的Scott-Zhang插值对 $H^1$ 中的函数有定义，将被用来推导有限元解的后验误差估计.

记 $\mathcal{N}_h$ 为网格 $\mathcal{M}_h$ 中节点集合. 对任意 $z \in \mathcal{N}_h$ , 记 $\phi_z \in V_h$ 为 $z$ 点对应的节点基函数. 显然,

$$
v _ {h} = \sum_ {z \in \mathcal {N} _ {h}} v _ {h} (z) \phi_ {z}, \quad \forall v _ {h} \in V _ {h}.
$$

取 $e_z$ 为以节点 $z$ 为顶点的一个单元的边/面, 要求其满足: 若 $z \in \partial \Omega$ , 则 $e_z \subset \partial \Omega$ .

引入 $e_z$ 上的线性函数：

$$
\psi_ {z} \in P _ {1} \left(e _ {z}\right): \int_ {e _ {z}} \psi_ {z} q d s = q (z), \quad \forall q \in P _ {1} \left(e _ {z}\right). \tag {11.4}
$$

显然

$$
\| \psi_ {z} \| _ {L ^ {2} (e _ {z})} ^ {2} = \int_ {e _ {z}} \psi_ {z} \psi_ {z} d s = \psi_ {z} (z) \leqslant C | e _ {z} | ^ {- \frac {1}{2}} \| \psi_ {z} \| _ {L ^ {2} (e _ {z})},
$$

从而 $\psi_z$ 满足估计

$$
\| \psi_ {z} \| _ {L ^ {2} \left(e _ {z}\right)} \leqslant C | e _ {z} | ^ {- \frac {1}{2}}, \quad \| \psi_ {z} \| _ {L ^ {\infty} \left(e _ {z}\right)} \leqslant C | e _ {z} | ^ {- 1}. \tag {11.5}
$$

任给 $v \in H^{1}(\Omega)$ , 定义其Scott-Zhang插值 $\Pi_h v \in V_h$ 满足:

$$
(\varPi_ {h} v) (z) = \int_ {e _ {z}} \psi_ {z} v \mathrm {d} s, \quad \forall z \in \mathcal {N} _ {h}, \quad \text {即} \quad \varPi_ {h} v = \sum_ {z \in \mathcal {N} _ {h}} \left(\int_ {e _ {z}} \psi_ {z} v \mathrm {d} s\right) \phi_ {z}. \qquad (1 1. 6)
$$

下面定理说明Scott-Zhang插值是一个投影，自然保持齐次Dirichlet边界条件，并给出其误差估计.

定理11.1 存在一个只依赖于 $\mathcal{M}_h$ 的最小角的常数 $C$ ，使得对于任意 $v\in H^{1}(\Omega)$ $K\in \mathcal{M}_h$ 有

(i) $\Pi_h v_h = v_h, \forall v_h \in V_h$ ，即 $\Pi_h$ 是一个投影算子.  
(ii) $\varPi_{h}v\in V_{h}^{0},\forall v\in H_{0}^{1}(\varOmega).$   
(iii) $\| v - \Pi_h v\|_{L^2 (K)} + h_K\| \nabla (v - \Pi_h v)\|_{L^2 (K)}\leqslant Ch_K\| \nabla v\|_{L^2 (\tilde{K})}.$   
(iv) $\| v - \Pi_h v \|_{L^2 (\partial K)} \leqslant Ch_e^{\frac{1}{2}} \| \nabla v \|_{L^2 (\tilde{K})}$ .   
(v) $\| \nabla \Pi_h v\|_{L^2 (K)}\leqslant C\| \nabla v\|_{L^2 (\tilde{K})},$

$$
\left\| \Pi_ {h} v \right\| _ {L ^ {2} (K)} \leqslant C \left(\left\| v \right\| _ {L ^ {2} (K)} + h _ {K} \left\| \nabla v \right\| _ {L ^ {2} (\tilde {K})}\right).
$$

其中 $\tilde{K}$ 是 $\mathcal{M}_h$ 中与 $K$ 的交非空的所有单元的并集

证明 由 $\varPi_h$ 的定义，(i)和(ii)显然成立．(iv)可以由(iii）和局部迹不等式(6.77)推出.(v)显然是(iii)和三角不等式的推论.剩下只需证明(iii).先证明如下稳定性估计：

$$
\left\| \Pi_ {h} v \right\| _ {L ^ {2} (K)} + h _ {K} \left\| \nabla \Pi_ {h} v \right\| _ {L ^ {2} (K)} \leqslant \left\| v \right\| _ {L ^ {2} (\tilde {K})} + h _ {K} \left\| \nabla v \right\| _ {L ^ {2} (\tilde {K})}. \tag {11.7}
$$

首先，对任一顶点 $z \in K$ ，记 $K_z$ 为以 $e_z$ 为其一边/面的单元，由（11.5）及局部迹不等式(6.77)得

$$
\left| \left(\Pi_ {h} v\right) (z) \right| \leqslant \| \psi_ {z} \| _ {L ^ {2} \left(e _ {z}\right)} \| v \| _ {L ^ {2} \left(e _ {z}\right)} \leqslant C | e _ {z} | ^ {- \frac {1}{2}} \left(h _ {K} ^ {- \frac {1}{2}} \| v \| _ {L ^ {2} \left(K _ {z}\right)} + h _ {K} ^ {\frac {1}{2}} \| \nabla v \| _ {L ^ {2} \left(K _ {z}\right)}\right).
$$

故由有限元逆不等式 (6.76) 得

$$
\begin{array}{l} \| \Pi_ {h} v \| _ {L ^ {2} (K)} + h _ {K} \| \nabla \Pi_ {h} v \| _ {L ^ {2} (K)} \leqslant C \| \Pi_ {h} v \| _ {L ^ {2} (K)} \leqslant C | K | ^ {\frac {1}{2}} \max  _ {z \in K} | (\Pi_ {h} v) (z) | \\ \leqslant C | K | ^ {\frac {1}{2}} \max  _ {z \in K} | e _ {z} | ^ {- \frac {1}{2}} \left(h _ {K} ^ {- \frac {1}{2}} \| v \| _ {L ^ {2} \left(K _ {z}\right)} + h _ {K} ^ {\frac {1}{2}} \| \nabla v \| _ {L ^ {2} \left(K _ {z}\right)}\right) \\ \leqslant C | K | ^ {\frac {1}{2}} | \partial K | ^ {- \frac {1}{2}} h _ {K} ^ {- \frac {1}{2}} \left(\| v \| _ {L ^ {2} (\bar {K})} + h _ {K} \| \nabla v \| _ {L ^ {2} (\bar {K})}\right), \\ \end{array}
$$

从而 (11.7) 成立. 由 (i) 和 (11.7) 得, 对任意 $v_{h} \in V_{h}$ 有

$$
\begin{array}{l} \| v - \Pi_ {h} v \| _ {L ^ {2} (K)} + h _ {K} \| \nabla (v - \Pi_ {h} v) \| _ {L ^ {2} (K)} \\ = \| v - v _ {h} - \Pi_ {h} (v - v _ {h}) \| _ {L ^ {2} (K)} + h _ {K} \| \nabla (v - v _ {h} - \Pi_ {h} (v - v _ {h})) \| _ {L ^ {2} (K)} \\ \leqslant C \left(\| v - v _ {h} \| _ {L ^ {2} (\tilde {K})} + h _ {K} \| \nabla (v - v _ {h}) \| _ {L ^ {2} (\tilde {K})}\right). \tag {11.8} \\ \end{array}
$$

注意到网格是正则的，以每个 $z \in \mathcal{N}_h$ 为顶点的单元的个数有上界. 将 $\tilde{K}$ 视为一个小网格片，那么存在有限个参考网格片 $\widehat{\tilde{K}}_1, \widehat{\tilde{K}}_2, \dots, \widehat{\tilde{K}}_m$ ，其中 $m$ 只与网格的最小角有关，使得对任一 $K \in \mathcal{M}_h$ ， $\tilde{K}$ 都能分片仿射等价于某个 $\widehat{\tilde{K}}_{i_K}, 1 \leqslant i_K \leqslant m.$ 由尺度变换技巧（见引理6.3—6.4），得

$$
\left\| v - v _ {h} \right\| _ {L ^ {2} (\tilde {K})} + h _ {K} \| \nabla (v - v _ {h}) \| _ {L ^ {2} (\tilde {K})} \leqslant C | K | ^ {\frac {1}{2}} \left(\left\| \hat {v} - \hat {v} _ {h} \right\| _ {L ^ {2} (\widehat {\tilde {K}})} + \left\| \widehat {\nabla} (\hat {v} - \hat {v} _ {h}) \right\| _ {L ^ {2} (\widehat {\tilde {K}})}\right).
$$

取 $v_{h}$ 使得 $\hat{v}_h$ 为 $\hat{v}$ 在 $\hat{\tilde{K}}$ 上的积分平均并利用Poincaré不等式得

$$
\left\| v - v _ {h} \right\| _ {L ^ {2} (\tilde {K})} + h _ {K} \left\| \nabla \left(v - v _ {h}\right) \right\| _ {L ^ {2} (\tilde {K})} \leqslant C | K | ^ {\frac {1}{2}} \left\| \widehat {\nabla} \hat {v} \right\| _ {L ^ {2} (\widehat {\tilde {K}})} \leqslant C h _ {K} \left\| \nabla v \right\| _ {L ^ {2} (\tilde {K})}.
$$

代入 (11.8) 即得 (iii) 成立. 证毕

# 11.2.2 后验误差估计

对于任意 $e \in \mathcal{E}_h^I$ 且 $e = K_1 \cap K_2$ , 我们定义 $u_h$ 的跳跃残量

$$
J _ {e} = \llbracket a (x) \nabla u _ {h} \rrbracket \big | _ {e} := a (x) \nabla u _ {h} \mid_ {K _ {1}} \cdot \nu_ {1} + a (x) \nabla u _ {h} \mid_ {K _ {2}} \cdot \nu_ {2} \tag {11.9}
$$

其中 $\nu_{i}$ 是 $\partial K_{i}$ 的单位外法向量在 $e$ 上的限制.为了方便，对任意一边界边/面 $e\subset \partial \Omega$ 定义 $J_{e} = 0.$ 对于任意单元 $K\in \mathcal{M}_h$ ，定义误差指示子 $\eta_K$ 为

$$
\eta_ {K} ^ {2} := h _ {K} ^ {2} \| f \| _ {L ^ {2} (K)} ^ {2} + h _ {K} \sum_ {e \subset \partial K} \| J _ {e} \| _ {L ^ {2} (e)} ^ {2}. \tag {11.10}
$$

显然，一旦算出了有限元解 $u_{h}$ ，那么 $\eta_K$ 是可计算的

对任意子网格 $\mathcal{T} \subset \mathcal{M}_h$ ，记 $\eta_{\mathcal{T}} = \left(\sum_{K \in \mathcal{T}} \eta_K^2\right)^{\frac{1}{2}}$ 。对任何子区域 $G \subset \Omega$ ，令 $\|\cdot\|_G = \|a^{\frac{1}{2}} \nabla \cdot\|_{L^2(G)}$ 。注意到， $\|\cdot\|_{\Omega}$ 即为 $H_0^1(\Omega)$ 上的能量范数。下面定理给出有限元解误差的一个后验的上界估计。

定理11.2 存在一个只依赖于网格 $\mathcal{M}_h$ 的最小角度和 $a(x)$ 最小值的常数 $C_1 > 0$ 使得

$$
\| u - u _ {h} \| _ {\Omega} \leqslant C _ {1} \eta_ {\mathcal {M} _ {h}}.
$$

证明 定义有限元解的残量 $R \in H^{-1}(\Omega)$

$$
\langle R, \varphi \rangle = (f, \varphi) - (a \nabla u _ {h}, \nabla \varphi) = (a \nabla (u - u _ {h}), \nabla \varphi), \quad \forall \varphi \in H _ {0} ^ {1} (\Omega).
$$

通过 (11.2) (11.3) 我们得到 Galerkin 正交性: $\langle R, u_h \rangle = 0$ , $\forall v_h \in V_h^0$ . 因此,

$$
\begin{array}{l} (a \nabla (u - u _ {h}), \nabla \varphi) = \langle R, \varphi - \Pi_ {h} \varphi \rangle \\ = (f, \varphi - \Pi_ {h} \varphi) - (a \nabla u _ {h}, \nabla (\varphi - \Pi_ {h} \varphi)) \\ = (f, \varphi - \Pi_ {h} \varphi) - \sum_ {K \in \mathcal {M} _ {h}} \int_ {K} a \nabla u _ {h} \cdot \nabla (\varphi - \Pi_ {h} \varphi) d x \\ = (f, \varphi - \Pi_ {h} \varphi) - \sum_ {K \in \mathcal {M} _ {h}} \int_ {\partial K} a \nabla u _ {h} \cdot \nu (\varphi - \Pi_ {h} \varphi) d x \\ = \sum_ {K \in \mathcal {M} _ {h}} \int_ {K} f (\varphi - \Pi_ {h} \varphi) d x - \sum_ {e \in \mathcal {E} _ {h} ^ {I}} \int_ {e} J _ {e} (\varphi - \Pi_ {h} \varphi) d x \\ \leqslant C \left(\sum_ {K \in \mathcal {M} _ {h}} h _ {K} ^ {2} \| f \| _ {L ^ {2} (K)} ^ {2}\right) ^ {\frac {1}{2}} \| \nabla \varphi \| _ {L ^ {2} (\Omega)} + \\ C \left(\sum_ {e \in \mathcal {E} _ {h} ^ {I}} h _ {e} \| J _ {e} \| _ {L ^ {2} (e)} ^ {2}\right) ^ {\frac {1}{2}} \| \nabla \varphi \| _ {L ^ {2} (\Omega)} \\ \leqslant C _ {1} \left(\sum_ {K \in \mathcal {M} _ {h}} \eta_ {K} ^ {2}\right) ^ {\frac {1}{2}} \| \vert \varphi \vert \| _ {\Omega}. \\ \end{array}
$$

取 $\varphi = u - u_{h}\in H_{0}^{1}(\Omega)$ 即得证明

下面定理给出局部的下界估计

定理11.3 存在一个只依赖于网格 $\mathcal{M}_h$ 的最小角度和 $a(x)$ 最大值的常数 $C_2 > 0$ 使得对任意单元 $K\in \mathcal{M}_h$ 有

$$
\eta_ {K} ^ {2} \leqslant C _ {2} \| u - u _ {h} \| _ {K ^ {*}} ^ {2} + C _ {2} \sum_ {K \subset K ^ {*}} h _ {K} ^ {2} \| f - f _ {K} \| _ {L ^ {2} (K)} ^ {2},
$$

其中 $f_{K} = \frac{1}{|K|}\int_{K}f\mathrm{d}\pmb {x},K^{*}$ 是与 $K$ 至少有一公共边/面的所有单元的并集

证明 由定理11.2的证明，

$$
(a \nabla (u - u _ {h}), \nabla \varphi) = \sum_ {K \in \mathcal {M} _ {h}} \int_ {K} f \varphi \mathrm {d} x - \sum_ {e \in \mathcal {E} _ {h} ^ {I}} \int_ {e} J _ {e} \varphi \mathrm {d} s, \quad \forall \varphi \in H _ {0} ^ {1} (\Omega). \tag {11.11}
$$

剩下的证明分为两步

$1^{\circ}$ 对任意 $K\in \mathcal{M}_h$ ，令 $\varphi_{K} = (d + 1)^{d + 1}\lambda_{1}\lambda_{2}\dots \lambda_{d + 1}$ 为 $K$ 上的标准泡函数，取 $\varphi = \left\{ \begin{array}{ll}\varphi_Kf_K, & \boldsymbol {x}\in K,\\ 0, & \boldsymbol {x}\notin K, \end{array} \right.$ 则易知

$$
\| f _ {K} \| _ {L ^ {2} (K)} ^ {2} \leqslant C \int_ {K} f _ {K} \varphi d \boldsymbol {x},
$$

$$
\| \varphi \| _ {L ^ {2} (K)} + h _ {K} \| \nabla \varphi \| _ {L ^ {2} (K)} \leqslant C \| \varphi \| _ {L ^ {2} (K)} \leqslant C \| f _ {K} \| _ {L ^ {2} (K)}.
$$

由 (11.11) 可得

$$
\begin{array}{l} \| f _ {K} \| _ {L ^ {2} (K)} ^ {2} \leqslant C \int_ {K} f _ {K} \varphi \mathrm {d} \boldsymbol {x} = C \left(\int_ {K} \left(f _ {K} - f\right) \varphi \mathrm {d} \boldsymbol {x} + \int_ {K} a \nabla \left(u - u _ {h}\right) \nabla \varphi \mathrm {d} \boldsymbol {x}\right) \\ \leqslant C \| f - f _ {K} \| _ {L ^ {2} (K)} \| \varphi \| _ {L ^ {2} (K)} + C \| | u - u _ {h} | \| _ {K} \| \nabla \varphi \| _ {L ^ {2} (K)} \\ \leqslant C \| f _ {K} \| _ {L ^ {2} (K)} \left(h _ {K} ^ {- 1} \| u - u _ {h} \| _ {K} + \| f - f _ {K} \| _ {L ^ {2} (K)}\right). \\ \end{array}
$$

因此,

$$
h _ {K} \| f _ {K} \| _ {L ^ {2} (K)} \leqslant C \left(\| u - u _ {h} \| _ {K} + h _ {K} \| f - f _ {K} \| _ {L ^ {2} (K)}\right).
$$

再由三角不等式得

$$
\left\| h _ {K} f \right\| _ {L ^ {2} (K)} ^ {2} \leqslant C \left(\left\| u - u _ {h} \right\| _ {K} ^ {2} + \left\| h _ {K} (f - f _ {K}) \right\| _ {L ^ {2} (K)} ^ {2}\right).
$$

$2^{\circ}$ 对于任意边/面 $e \subset \partial K \cap \Omega$ ，记 $\Omega_{e}$ 为以 $e$ 为公共边/面的两个单元的并，设 $\psi_{e} = d^{d}\lambda_{1}\lambda_{2}\dots \lambda_{d}$ 为 $\Omega_{e}$ 上的泡函数，其中 $\lambda_1,\lambda_2,\dots ,\lambda_d$ 是与 $e$ 的节点对应的重心坐标函数.取 $\psi = \left\{ \begin{array}{ll}\psi_{e}J_{e}, & x\in \Omega_{e},\\ 0, & x\notin \Omega_{e}, \end{array} \right.$ 则

$$
\| J _ {e} \| _ {L ^ {2} (e)} ^ {2} \leqslant C \int_ {e} J _ {e} \psi \mathrm {d} s,
$$

容易验证因此

$$
\| \psi \| _ {L ^ {2} \left(\Omega_ {e}\right)} + h _ {K} \| \nabla \psi \| _ {L ^ {2} \left(\Omega_ {e}\right)} \leqslant C \| \psi \| _ {L ^ {2} \left(\Omega_ {e}\right)} \leqslant C h _ {e} ^ {\frac {1}{2}} \| J _ {e} \| _ {L ^ {2} (e)}.
$$

现在由 (11.11) 和 $\psi \in H_0^1 (\Omega_e)$ 得

$$
\begin{array}{l} \| J _ {e} \| _ {L ^ {2} (e)} ^ {2} \leqslant C \int_ {e} J _ {e} \psi \mathrm {d} s = C \left(\int_ {\Omega_ {e}} f \psi \mathrm {d} x - \int_ {\Omega_ {e}} a \nabla (u - u _ {h}) \nabla \psi \mathrm {d} x\right) \\ \leqslant C \| J _ {e} \| _ {L ^ {2} (e)} \left(h _ {e} ^ {\frac {1}{2}} \| f \| _ {L ^ {2} (\Omega_ {e})} + h _ {e} ^ {- \frac {1}{2}} \| u - u _ {h} \| _ {\Omega_ {e}}\right). \\ \end{array}
$$

所以代入 $\| h_K f\|_{L^2 (K)}$ 的估计即得证明

定理11.3表明，若 $f$ 分片光滑，则在相差一个高阶量 $\left(\sum_{K\subset K^{\star}}h_K^2\| f - f_K\|_{L^2 (K)}^2\right)^{\frac{1}{2}}$ 的前提下，误差指示子 $\eta_{K}$ 可以作为局部能量误差 $\| u - u_h\|_{K^*}$ 的下界.所以，如果 $\eta_{K}$ 大，那么在单元 $K$ 附近误差一定大.

# 11.3 自适应算法

基于局部误差指示子，求解变分问题（11.3）的自适应算法通常可以描述为如下形式的循环：

$$
\text {求 解} \longrightarrow \text {估 计} \longrightarrow \text {标 记} \longrightarrow \text {加 密}. \tag {11.12}
$$

为了保证自适应算法收敛，即，从任何给定初始网格开始，循环迭代（11.12）有限步终止，我们需要适当设计标记策略。文献中已有的标记策略，通常基于所谓的误差平均分配原则：最优的网格应该按单元平均分配误差。下面，我们简要回顾两种常用的标记策略：已知网格 $\mathcal{M}_H$ 和其上的误差指示子 $\eta_K, K \in \mathcal{M}_H$

1. 最大值策略：给定 $\theta \in (0,1)$ ，标记所有满足

$$
\eta_ {K} \geqslant \theta \max  _ {K ^ {\prime} \in \mathcal {M} _ {H}} \eta_ {K ^ {\prime}}
$$

的单元 $K$

2. Dörfler 策略: 给定 $\theta \in (0,1]$ , 找 $\mathcal{M}_H$ 的子集 $\hat{\mathcal{M}}_H$ 并标记其中的单元, 使得

$$
\eta_ {\mathcal {M} _ {H}} \geqslant \theta \eta_ {\mathcal {M} _ {H}}. \tag {11.13}
$$

需要说明的是，在Dörfler策略中，为了保证自适应算法的拟最优性，我们一般选取最少个数的单元使得（11.13）成立。

给定一个粗网格 $\mathcal{M}_H$ 和标记的子网格 $\hat{\mathcal{M}}_H\subset \mathcal{M}_H$ ，通过加密 $\hat{\mathcal{M}}_H$ 中的单元得到细网格 $\mathcal{M}_h$ .当然，加密网格通常包括两个步骤：加密标记的单元和去悬点（hangingnodes).我们对第一步作如下假设：存在常数 $m > 1$ ，使得

$$
| K | \leqslant \frac {1}{m} | K ^ {\prime} |, \quad \forall K \subset K ^ {\prime}, K \in \mathcal {M} _ {h}, K ^ {\prime} \in \hat {\mathcal {M}} _ {H}, \tag {11.14}
$$

即任何标记的粗网格单元加密后所得细网格子单元的测度不超过这个粗网格单元测度的 $\frac{1}{m}$ . 例如, 在二等分加密的情况下, $m = 2$ . 注意到, 在去悬点的步骤中, 可能会加密一些未标记的单元.

下面给出自适应有限元算法11.1：

# 算法11.1 自适应有限元算法

给定 $\theta \in (0,1]$ ，初始网格 $\mathcal{M}_0$ ，赋值 $k\gets 0$

1. 求 $\mathcal{M}_k$ 上的有限元解 $u_{k}$   
2. 计算 $\mathcal{M}_k$ 上误差指示子 $\eta_K, \forall K \in \mathcal{M}_k$ .  
3. 按 Dörfler 策略标记单元, 即选 $\hat{\mathcal{M}}_k \subset \mathcal{M}_k$ 使得 $\eta_{\hat{\mathcal{M}}_k} \geqslant \theta \eta_{\mathcal{M}_k}$ .   
4. 加密 $\hat{\mathcal{M}}_k$ 得到 $\mathcal{M}_k$ , 满足 (11.14).  
5. $k \gets k + 1$ , 回到步 1.

# 11.4 收敛性分析

在本节中, 我们考虑基于 Dörfler 策略的自适应有限元算法 11.1 的收敛性. 我们从下面的引理开始.

引理11.2 设 $\mathcal{M}_h$ 是 $\mathcal{M}_H$ 的加密使得 $V_{H}\subset V_{h}$ ，则

$$
\left\| \left| u - u _ {h} \right| \right\| _ {\Omega} ^ {2} = \left\| \left| u - u _ {H} \right| \right\| _ {\Omega} ^ {2} - \left\| \left| u _ {h} - u _ {H} \right| \right\| _ {\Omega} ^ {2}.
$$

证明 由 $u_{h} - u_{H}\in V_{h}^{0}$ 及Galerkin正交性易证

令 $\widetilde{h}_K\coloneqq |K|^{\frac{1}{d}}$ ，引入与 $\eta_{K}$ 等价的误差指示子：

$$
\widetilde {\eta} _ {K} ^ {2} := \widetilde {h} _ {K} ^ {2} \| f \| _ {L ^ {2} (K)} ^ {2} + \widetilde {h} _ {K} \sum_ {e \subset \partial K} \| J _ {e} \| _ {L ^ {2} (e)} ^ {2}. \tag {11.15}
$$

显然，存在正常数 $c_{1}$ 和 $c_{2}$ ，使得

$$
c _ {2} \eta_ {K} \leqslant \widetilde {\eta} _ {K} \leqslant c _ {1} \eta_ {K}. \tag {11.16}
$$

修改后的误差指标 $\widetilde{\eta}_K$ 具有以下的缩减性质

引理11.3 设 $\hat{\mathcal{M}}_H \subset \mathcal{M}_H$ 为标记的粗网格单元的集合, $\mathcal{M}_h$ 为满足假设(11.14)的 $\mathcal{M}_H$ 的一个加密. 那么存在一个只依赖于网格 $\mathcal{M}_h$ 的最小角度和 $a(x)$ 最大值的一个常数 $C_3$ , 使得对任何 $\delta > 0$ , 有

$$
\widetilde {\eta} _ {\mathcal {M} _ {h}} ^ {2} \leqslant (1 + \delta) \left(\widetilde {\eta} _ {\mathcal {M} _ {H}} ^ {2} - \left(1 - \frac {1}{\sqrt [ d ]{m}}\right) \widetilde {\eta} _ {\mathcal {M} _ {H}} ^ {2}\right) + \left(1 + \frac {1}{\delta}\right) C _ {3} \| u _ {h} - u _ {H} \| _ {\Omega} ^ {2}.
$$

证明 由参数取为 $\delta$ 的 Young 不等式,

$$
\begin{array}{l} \widetilde {\eta} _ {\mathcal {M} _ {h}} ^ {2} = \sum_ {K \in \mathcal {M} _ {h}} \left(\widetilde {h} _ {K} ^ {2} \| f \| _ {L ^ {2} (K)} ^ {2} + \widetilde {h} _ {K} \sum_ {e \subset \partial K \cap \Omega} \| [ a \nabla (u _ {H} + u _ {h} - u _ {H}) ] \| _ {L ^ {2} (e)} ^ {2}\right) \\ \leqslant \sum_ {K \in \mathcal {M} _ {h}} \left(\widetilde {h} _ {K} ^ {2} \| f \| _ {L ^ {2} (K)} ^ {2} + (1 + \delta) \widetilde {h} _ {K} \sum_ {e \subset \partial K \cap \Omega} \| [ a \nabla u _ {H} ] \| _ {L ^ {2} (e)} ^ {2}\right) + \\ \left(1 + \frac {1}{\delta}\right) \sum_ {K \in \mathcal {M} _ {h}} \widetilde {h} _ {K} \sum_ {e \subset \partial K \cap \Omega} \| [ a \nabla (u _ {h} - u _ {H}) ] \| _ {L ^ {2} (e)} ^ {2} \\ =: I + I I. \\ \end{array}
$$

注意到, 对于一个粗网格单元 $K' \in \mathcal{M}_H$ 内部的任何 $e$ 有 $[a\nabla u_H]|_e = 0$ , 并且对任何细网格单元 $K \subset K' \in \hat{\mathcal{M}}_H$ 有 $\widetilde{h}_K = |K|^{\frac{1}{d}} \leqslant \frac{1}{\sqrt[d]{m}} \widetilde{H}_{K'}$ , 我们有

$$
I \leqslant (1 + \delta) \sum_ {K \subset K ^ {\prime} \in \mathcal {M} _ {H} \backslash \hat {\mathcal {M}} _ {H}} \left(\widetilde {h} _ {K} ^ {2} \| f \| _ {L ^ {2} (K)} ^ {2} + \widetilde {h} _ {K} \sum_ {e \subset \partial K \cap \Omega} \| [ a \nabla u _ {H} ] \| _ {L ^ {2} (e)} ^ {2}\right) +
$$

$$
\begin{array}{l} (1 + \delta) \sum_ {K \subset K ^ {\prime} \in \hat {\mathcal {M}} _ {H}} \left(\widetilde {h} _ {K} ^ {2} \| f \| _ {L ^ {2} (K)} ^ {2} + \widetilde {h} _ {K} \sum_ {e \subset \partial K \cap \Omega} \| [ a \nabla u _ {H} ] \| _ {L ^ {2} (e)} ^ {2}\right) \\ \leqslant (1 + \delta) \sum_ {K ^ {\prime} \in \mathcal {M} _ {H} \backslash \hat {\mathcal {M}} _ {H}} \left(\widetilde {H} _ {K ^ {\prime}} ^ {2} \| f \| _ {L ^ {2} \left(K ^ {\prime}\right)} ^ {2} + \widetilde {H} _ {K ^ {\prime}} \sum_ {e ^ {\prime} \subset \partial K ^ {\prime} \cap \Omega} \| [ a \nabla u _ {H} ] \| _ {L ^ {2} \left(e ^ {\prime}\right)} ^ {2}\right) + \\ \frac {1 + \delta}{\sqrt [ d ]{m}} \sum_ {K ^ {\prime} \in \tilde {\mathcal {M}} _ {H}} \left(\widetilde {H} _ {K ^ {\prime}} ^ {2} \| f \| _ {L ^ {2} \left(K ^ {\prime}\right)} ^ {2} + \widetilde {H} _ {K ^ {\prime}} \sum_ {e ^ {\prime} \subset \partial K ^ {\prime} \cap \Omega} \| [ a \nabla u _ {H} ] \| _ {L ^ {2} \left(e ^ {\prime}\right)} ^ {2}\right) \\ = (1 + \delta) \widetilde {\eta} _ {\mathcal {M} _ {H} \backslash \hat {\mathcal {M}} _ {H}} ^ {2} + \frac {1 + \delta}{\sqrt [ d ]{m}} \widetilde {\eta} _ {\hat {\mathcal {M}} _ {H}} ^ {2} = (1 + \delta) \left(\widetilde {\eta} _ {\mathcal {M} _ {H}} ^ {2} - \left(1 - \frac {1}{\sqrt [ d ]{m}}\right) \widetilde {\eta} _ {\hat {\mathcal {M}} _ {H}} ^ {2}\right). \\ \end{array}
$$

接下来我们估计II.对任意 $e\in \mathcal{E}_h^I$ ，记 $K_{1}^{e}$ 和 $K_{2}^{e}$ 为以 $e$ 为公共边/面的两个单元我们有

$$
\begin{array}{l} I I \leqslant C \left(1 + \frac {1}{\delta}\right) \sum_ {e \in \mathcal {E} _ {h} ^ {I}} h _ {e} \| [ a \nabla (u _ {h} - u _ {H}) ] \| _ {L ^ {2} (e)} ^ {2} \\ \leqslant C \left(1 + \frac {1}{\delta}\right) \sum_ {e \in \mathcal {E} _ {h} ^ {I}} h _ {e} \left(\| a \nabla \left(u _ {h} - u _ {H}\right) | _ {K _ {1} ^ {e}} \| _ {L ^ {2} (e)} ^ {2} + \| a \nabla \left(u _ {h} - u _ {H}\right) | _ {K _ {2} ^ {e}} \| _ {L ^ {2} (e)} ^ {2}\right) \\ \leqslant C \left(1 + \frac {1}{\delta}\right) \sum_ {K \in \mathcal {M} _ {h}} h _ {K} \| a \nabla (u _ {h} - u _ {H}) \| _ {L ^ {2} (\partial K)} ^ {2} \\ \leqslant \left(1 + \frac {1}{\delta}\right) C _ {3} \| \left| u _ {h} - u _ {H} \right| \| _ {\Omega} ^ {2}, \\ \end{array}
$$

这里用到了局部迹不等式 (6.77) 来推导最后一个不等式. 结合上面三个估计即得证明

下面定理说明，自适应有限元迭代得到的有限元解序列的误差按线性速度趋于零。

定理11.4 给定 $\theta \in (0,1]$ , 设 $\{\mathcal{M}_k, u_k\}_{k \geqslant 0}$ 为按自适应有限元算法11.1产生的网格和离散解序列. 假设网格族 $\{\mathcal{M}_k\}$ 是正则的, 则存在常数 $\gamma > 0$ , $C_0 > 0$ 和 $0 < \alpha < 1$ , 仅依赖于 $\{\mathcal{M}_k\}$ 的形状正则性, $m$ 和标记参数 $\theta$ , 使得

$$
\left(\| \| u - u _ {k} \| \| _ {\Omega} ^ {2} + \gamma \eta_ {\mathcal {M} _ {k}} ^ {2}\right) ^ {\frac {1}{2}} \leqslant C _ {0} \alpha^ {k}. \tag {11.17}
$$

证明 我们首先证明存在常数 $\gamma_0 > 0$ 和 $0 < \alpha < 1$ , 使得

$$
\left\| \left\| u - u _ {k + 1} \right\| \right\| _ {\Omega} ^ {2} + \gamma_ {0} \widetilde {\eta} _ {\mathcal {M} _ {k + 1}} ^ {2} \leqslant \alpha^ {2} \left(\left\| \left\| u - u _ {k} \right\| \right\| _ {\Omega} ^ {2} + \gamma_ {0} \widetilde {\eta} _ {\mathcal {M} _ {k}} ^ {2}\right). \tag {11.18}
$$

为了简便，记

$$
e _ {k} := \| \left| u - u _ {k} \right| \| _ {\Omega}, \quad \widetilde {\eta} _ {k} := \widetilde {\eta} _ {\mathcal {M} _ {k}}, \quad \lambda := 1 - \frac {1}{\sqrt [ d ]{m}}.
$$

由引理11.2—11.3和Dörfler策略，我们有

$$
\widetilde {\eta} _ {k + 1} \leqslant (1 + \delta) (1 - \lambda \theta^ {2}) \widetilde {\eta} _ {k} + \left(1 + \frac {1}{\delta}\right) C _ {3} \left(e _ {k} ^ {2} - e _ {k + 1} ^ {2}\right). \tag {11.19}
$$

由上界估计定理11.2和(11.16)，我们有

$$
e _ {k} ^ {2} \leqslant \widetilde {C} _ {1} \widetilde {\eta} _ {k} ^ {2}, \quad \text {其 中} \quad \widetilde {C} _ {1} = \frac {C _ {1}}{c _ {2}}. \tag {11.20}
$$

记 $\beta_{\mathrm{i}} = \left(1 + \frac{1}{\delta}\right)C_{3}.$ 对任意 $0 <   \zeta <  1$ ，由（11.19）一(11.20）知，

$$
\begin{array}{l} e _ {k + 1} ^ {2} + \frac {1}{\beta} \widetilde {\eta} _ {k + 1} ^ {2} \leqslant e _ {k} ^ {2} + \frac {1}{\beta} (1 + \delta) (1 - \lambda \theta^ {2}) \widetilde {\eta} _ {k} ^ {2} \\ \leqslant \zeta e _ {k} ^ {2} + \left[ (1 - \zeta) \widetilde {C} _ {1} + \frac {1}{\beta} (1 + \delta) (1 - \lambda \theta^ {2}) \right] \widetilde {\eta} _ {k} ^ {2} \\ = \zeta \left[ e _ {k} ^ {2} + \frac {1}{\beta} \left(\beta \zeta^ {- 1} (1 - \zeta) \widetilde {C} _ {1} + \zeta^ {- 1} (1 + \delta) (1 - \lambda \theta^ {2})\right) \widetilde {\eta} _ {k} ^ {2} \right]. \\ \end{array}
$$

我们取 $\delta > 0$ 使得 $(1 + \delta)(1 - \lambda \theta^2) < 1$ , 并取 $\zeta$ 使得 $\beta \zeta^{-1}(1 - \zeta)\widetilde{C}_1 + \zeta^{-1}(1 + \delta)(1 - \lambda \theta^2) = 1$ , 即取

$$
\zeta = \frac {(1 + \delta) (1 - \lambda \theta^ {2}) + \beta \widetilde {C} _ {1}}{1 + \beta \widetilde {C} _ {1}} <   1,
$$

从而得 (11.18) 成立, 其中

$$
\gamma_ {0} = \frac {1}{\beta} = \frac {\delta}{(1 + \delta) C _ {3}}, \quad \alpha^ {2} = \zeta = \frac {\delta (1 + \delta) (1 - \lambda \theta^ {2}) + (1 + \delta) \widetilde {C} _ {1} C _ {3}}{\delta + (1 + \delta) \widetilde {C} _ {1} C _ {3}}.
$$

进一步, 注意到, 由 (11.16), $\widetilde{\eta}_k \geqslant c_2 \eta_k$ , 令

$$
\gamma = \gamma_ {0} c _ {2} ^ {2} \quad \text {和} \quad C _ {0} = \left(\| \| u - u _ {0} \| \| _ {\Omega} ^ {2} + \gamma_ {0} \widetilde {\eta} _ {0} ^ {2}\right) ^ {\frac {1}{2}}.
$$

即知 (11.17) 成立. 得证

在二维情况下, 大量的数值实验表明, 本章中描述的基于后验误差估计的自适应有限元方法不仅收敛, 而且得到网格和相关的数值复杂性是拟最优的, 即其上的线性有限元离散的能量范数误差为 $O(N^{-\frac{1}{2}})$ , 其中 $N$ 是自由度数. 这一点在适当条件下可以给出证明, 但由于篇幅所限, 我们略去自适应有限元算法的拟最优性理论.

例11.2 继续考虑例11.1中的 $L$ 形区域问题，这次我们用线性自适应有限元方法11.1求解。图11.3绘制了10次自适应迭代后的网格图(a)，和对数-对数坐标下的误差 $\| u - u_{k}\|_{H^{1}(\Omega)}$ 关于 $N_{k}$ 的曲线(b)，其中 $u_{k}$ 是 $k$ 次迭代后的网格 $\mathcal{M}_k$ 上的有限元解， $\mathcal{N}_k$ 是 $\mathcal{M}_k$ 上的总自由度数。可以看出

$$
\left\| u - u _ {k} \right\| _ {H ^ {1} (\Omega)} = O \left(N _ {k} ^ {- \frac {1}{2}}\right), \tag {11.21}
$$

即收敛速度是拟最优的

![](images/caa6495e4759481382518a2a2f494de3548149de00f6cdefaef974e9eda3c23e.jpg)  
(a) 经过10次自适应迭代后的网格

![](images/b981498b3fd3e2f4fa267cfef15ef464e76f7056fce4e81eff52f8672339dfe4.jpg)  
(b) 对数-对数坐标下的 $H^{1}$ 误差关于总自由度数的曲线. 虚线为斜率等于 $-\frac{1}{2}$ 的参考线  
图11.3

# 11.4.1 习题

1. 求扇形区域 $S_w$ 上的 Laplace 方程 $-\Delta u = 0$ 的形如 $u = r^{\alpha} \mu(\theta)$ 的通解, 分别满足边界条件:

(i) $\left.\frac{\partial u}{\partial v}\right|_{\Gamma_i} = 0, i = 1,2;$   
(ii) $u|_{\Gamma_1} = 0, \left.\frac{\partial u}{\partial v}\right|_{\Gamma_2} = 0.$

2. 设 $\Omega$ 为 $\mathbb{R}^d (d = 2,3)$ 中的有界多面体域. 证明以下 Scott-Zhang 插值算子的误差估计.

$$
\left\| \varphi - \Pi_ {h} \varphi \right\| _ {H ^ {k} (K)} \leqslant C h _ {K} ^ {2 - k} | \varphi | _ {H ^ {2} (\tilde {K})}, \quad \forall \varphi \in H ^ {2} (\Omega), \quad k = 0, 1.
$$

3. 设 $\Omega \subset \mathbb{R}^2$ 是一个有界多边形区域. $f \in L^{2}(\Omega)$ 和 $g \in L^{2}(\partial \Omega)$ , 设 $u \in H^{1}(\Omega)$ 为如下椭圆型问题的弱解:

$$
- \Delta u = f, x \in \Omega ; \quad \frac {\partial u}{\partial n} + u = g, x \in \partial \Omega .
$$

推导其线性有限元解的后验误差上界估计

4. 对于两点边值问题: $-u'' = f, x \in (0,1)$ ; $u(0) = \alpha$ , $u'(1) = \beta$ , 推导其线性有限元解的后验误差上界和下界估计.

# 参考文献

[1] 李荣华, 刘播. 微分方程数值解法. 4版. 北京: 高等教育出版社, 2009.  
[2] 伍卓群, 李勇. 常微分方程. 2版. 北京: 高等教育出版社, 2023.  
[3] 江泽坚, 孙善利. 泛函分析. 2 版. 北京: 高等教育出版社, 2005.  
[4] 复旦大学数学系. 数学物理方程. 北京: 高等教育出版社, 1979.  
[5] 冯果忱, 黄明游. 数值分析: 上册. 北京: 高等教育出版社, 2007.  
[6] 黄明游, 冯果忱. 数值分析: 下册. 北京: 高等教育出版社, 2007.  
[7] 徐绪海, 朱方生. 刚性微分方程的数值方法. 武汉: 武汉大学出版社, 1997.  
[8] 黄明游. 发展方程的有限元方法. 上海: 上海科学技术出版社, 1988.  
[9] 李荣华. 偏微分方程数值解法. 2版. 北京: 高等教育出版社, 2010.  
[10] 李荣华. 边值问题的 Galerkin 有限元法. 北京: 科学出版社, 2005.  
[11] 胡建伟, 汤怀民. 微分方程的数值方法. 北京: 科学出版社, 1999.  
[12] 郭本瑜. 偏微分方程的差分法. 北京: 科学出版社, 1988.  
[13] 向新民. 谱方法的数值分析. 北京: 科学出版社, 2000.  
[14] 马驷良, 李荣华. 关于矩阵族一致有界的代数准则 (综合报告). 吉林大学自然科学学报, 1986 (1): 21-36.  
[15] 矢崎信男, 野木达夫. 发展方程的数值分析. 王宝兴, 殷广济, 雷光耀, 译. 北京: 人民教育出版社, 1982.  
[16] Samarskii, A.A. and Andreev, V.B.. 椭圆型方程数值方法. 武汉大学计算数学教研室, 译. 北京: 科学出版社, 1984.  
[17] Hackbusch, W.. 多重网格法. 林群, 等译. 北京: 科学出版社, 1998.  
[18] Chen, Z. and Wu, H.. Selected Topics in Finite Element Methods. Beijing: Beijing Science Press, 2010.   
[19] Li, R., Chen, Z., and Wu, W. Generalized Difference Methods for Differential Equations - Numerical Analysis of Finite Volume Methods. New York: Marcel Dekker, Inc., 2000.   
[20] Li, Y. and Li, R.. Generalized difference methods on arbitrary quadrilateral networks. J. Comp. Math., 1999, 17, 653-672.   
[21] Lv, J. and Li, Y., L² Error Estimate of The Finite Volume Element Methods on Quadrilateral Meshes. Adv. Comp. Math. 2010, 33, 129-148.   
[22] Bramble, J.H., Multigrid Methods. New York: Longman Scientific & Technical, 1993.   
[23] Courant, R. and Hilbert, D.. Methods of Mathematical Physics. Vol II. Hoboken: J. Wiley & Sons Inc., 1962.

[24] Ciarlet, D.G.. The Finite Element Method for Elliptic Problems. Amsterdam: North Holland, 1978.   
[25] Gear, C.W., Numerical Initial Value Problem in Ordinary Differential Equation. London: Prentice Hall, 1971.   
[26] Henrici, P., Discrete Variable Methods in Ordinary Differential Equations. Hoboken: J. Wiley & Sons Inc., 1962.   
[27] Leunbert, J.D.. Computational Methods in Ordinary Differential Equations. Hoboken: J. Wiley & Sons Inc., 1973.   
[28] Richtmyer, R.D. and Morton, K.W.. Difference Methods for Initial Value Problems. 2nd edition. Hoboken: J. Wiley & Sons Inc., 1967.   
[29] Strang, G. and Fix, G.J., An Analysis of The Finite Element Method. London: Prentice Hall, 1973.   
[30] Saul'yev, V.K.. Difference Methods for Solving Equation of Parabolic Type. Moscow: Moscow Press, 1960.   
[31] Thomas, J. W.. Numerical Partial Differential Equations. New York: Springer Verlag, 1999.   
[32] Wang, J. and Ye, X.. A Weak Galerkin Finite Element Method for Second Order Elliptic Problems. J. Comput. Appl. Math., 2013, 241, 103-115.   
[33] Mu, L., Wang, J., and Ye, X.. A Weak Galerkin Finite Element Method with Polynomial Reduction. J. Comp. Appl. Math., 2015, 285, 45-58.   
[34] X. Ye and S. Zhang. A Stabilizer-Free Weak Galerkin Finite Element Method on Polytopal Meshes. J. Comp. Appl. Math., 2020, 371, 112699, 9.   
[35] Brenner, S. and Scott, R.. The Mathematical Theory of Finite Element Methods. 3rd edition. New York: Springer, 2008.   
[36] Arnold, D.N., Brezzi, F., Cockburn, B., and Marini, L. D.. Unified Analysis of Discontinuous Galerkin Methods for Elliptic Problems. SIAM J. Numer. Anal. 2002, 39, 1749-1779.   
[37] Cockburn, B., Gopalakrishnan, J., and Sayas, F.. A Projection-Based Error Analysis of HDG Methods. Math. Comp. 2010, 79, 1351-1367.   
[38] Cockburn, B. and Shu, C.-W.. The Local Discontinuous Galerkin Method for Convection-Diffusion Systems. SIAM J. Numer. Anal. 1998, 35, 2440-2463.   
[39] Rivière, B.. Discontinuous Galerkin Methods for Solving Elliptic and Parabolic Equations: Theory and Implementation. New York: SIAM, 2008.

# 数学家简介

![](images/563717a3ba3f1a38452f047ed7cf851e97e13a3b2282aef6cabc883d6bbbb4ab.jpg)

# 郑重声明

高等教育出版社依法对本书享有专有出版权。任何未经许可的复制、销售行为均违反《中华人民共和国著作权法》，其行为人将承担相应的民事责任和行政责任；构成犯罪的，将被依法追究刑事责任。为了维护市场秩序，保护读者的合法权益，避免读者误用盗版书造成不良后果，我社将配合行政执法部门和司法机关对违法犯罪的单位和个人进行严厉打击。社会各界人士如发现上述侵权行为，希望及时举报，我社将奖励举报有功人员。

反盗版举报电话 (010) 58581999 58582371

反盗版举报邮箱 dd@hep.com.cn

通信地址 北京市西城区德外大街4号

高等教育出版社知识产权与法律事务部

邮政编码 100120

# 读者意见反馈

为收集对教材的意见建议，进一步完善教材编写并做好服务工作，读者可将对本教材的意见建议通过如下渠道反馈至我社。

咨询电话 400-810-0598

反馈邮箱 hepsci@pub.hep.cn

通信地址 北京市朝阳区惠新东街4号富盛大厦1座高等教育出版社理科事业部

邮政编码 100029

# 防伪查询说明

用户购书后刮开封底防伪涂层，使用手机微信等软件扫描二维码，会跳转至防伪查询网页，获得所购图书详细信息。

防伪客服电话 (010) 58582300

# 图书在版编目（CIP）数据

微分方程数值解法/李荣华，李永海，武海军编著北京：高等教育出版社，2025.4.--ISBN978-7-04-063698-7

I.0241.8

中国国家版本馆CIP数据核字第2025M1S558号

Weifen Fangcheng Shuzhi Jiefa

策划编辑 宋玉文

责任编辑 宋玉文

封面设计 王洋

版式设计 童丹

责任绘图 杨伟露

责任校对吕红颖

责任印制 赵义民

出版发行 高等教育出版社

社址北京市西城区德外大街4号

邮政编码 100120

购书热线 010-58581118

咨询电话 400-810-0598

网址http://www.hep.edu.cn

http://www.hep.com.cn

网上订购 http://www.hepmall.com.cn

http://www.hepmall.com

http://www.hepmall.cn

印刷北京盛通印刷股份有限公司

开本 $787\mathrm{mm}\times 1092\mathrm{mm}$ 1/16

印张 20.75

字数 390千字

版次2025年4月第1版

印次 2025年4月第1次印刷

定价 54.00元

本书如有缺页、倒页、脱页等质量问题，请到所购图书销售部门联系调换

版权所有 侵权必究

物料号 63698-00