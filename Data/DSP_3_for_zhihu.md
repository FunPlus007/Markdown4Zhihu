# 数字信号处理 从DTFT到FFT （三）

## 前言

在前面两篇文章，我们已经简单介绍了离散时域傅里叶变换和Z变换，但是，这两种对于离散信号的变换方式，得到的频率信号还是在时间上连续的，但是根据计算机存储处理的要求，我们希望能得到一组离散的数据，下面就来看看今天的离散傅里叶变换（DFT）。N点DFT公式如下。 <img src="https://www.zhihu.com/equation?tex=e^{-j\frac{2\pi}{N}}(W_n)" alt="e^{-j\frac{2\pi}{N}}(W_n)" class="ee_img tr_noresize" eeimg="1"> 为旋转因子

<img src="https://www.zhihu.com/equation?tex=X(k) = DFT[x(n)]= \sum_{n=0}^{N-1}x(n)e^{-j\frac{2\pi}{N}kn},k=0,1,...,N-1\\
x(n) = IDFT[X(k)] = \frac{1}{N}\sum_{n=0}^{N-1}X(k)e^{j\frac{2\pi}{N}kn},n=0,1,...,N-1
" alt="X(k) = DFT[x(n)]= \sum_{n=0}^{N-1}x(n)e^{-j\frac{2\pi}{N}kn},k=0,1,...,N-1\\
x(n) = IDFT[X(k)] = \frac{1}{N}\sum_{n=0}^{N-1}X(k)e^{j\frac{2\pi}{N}kn},n=0,1,...,N-1
" class="ee_img tr_noresize" eeimg="1">

## 1、DFT的物理意义

先来回顾一下DTFT，可以看到对于x(n)进行DTFT之后，得到的X在频率上还是连续的，在归一化的虚轴 <img src="https://www.zhihu.com/equation?tex=\omega" alt="\omega" class="ee_img tr_noresize" eeimg="1">  上以 <img src="https://www.zhihu.com/equation?tex=2\pi" alt="2\pi" class="ee_img tr_noresize" eeimg="1"> 为周期。这样显然是不利于计算机的处理的，因此引入DFT，对DTFT之后的信号在进行离散化采样。

<img src="https://www.zhihu.com/equation?tex=X(e^{jw}) = DTFT[x(n)] = \sum_{n=-\infty}^{+\infty} x(n)e^{-jwn} \\
" alt="X(e^{jw}) = DTFT[x(n)] = \sum_{n=-\infty}^{+\infty} x(n)e^{-jwn} \\
" class="ee_img tr_noresize" eeimg="1">


<img src="https://www.zhihu.com/equation?tex=X(z) = ZT[x(n)] = \sum_{n=-\infty}^{+\infty}x(n)z^{-n}\\
" alt="X(z) = ZT[x(n)] = \sum_{n=-\infty}^{+\infty}x(n)z^{-n}\\
" class="ee_img tr_noresize" eeimg="1">

DFT和DTFT、ZT之间的关系，**DFT是DTFT在 <img src="https://www.zhihu.com/equation?tex=[0,2\pi]" alt="[0,2\pi]" class="ee_img tr_noresize" eeimg="1"> 范围上的虚轴上的N点等间隔采样，是ZT在 <img src="https://www.zhihu.com/equation?tex=e^{jw}" alt="e^{jw}" class="ee_img tr_noresize" eeimg="1"> 的单位圆上的N点等间隔采样。** 关系式如下：

<img src="https://www.zhihu.com/equation?tex=X(k) = X(z)|_{z=e^{j\frac{2\pi}{N}k}}\\
X(k) = X(e^{jw})|_{w=\frac{2\pi}{N}k}
" alt="X(k) = X(z)|_{z=e^{j\frac{2\pi}{N}k}}\\
X(k) = X(e^{jw})|_{w=\frac{2\pi}{N}k}
" class="ee_img tr_noresize" eeimg="1">
所以DFT其实就是有限长序列傅里叶变换的有限点离散采样，从而实现频域离散化，目的就是为了方便计算机的处理。

关于DFT的性质，需要需要注意一个循环卷积，两个信号的时域的N点循环卷积即为频域的两个信号的N点DFT的乘积，FFT的优化计算是基于DFT的共轭对称，关于DFT的性质不再过多介绍。



## 2、频域采样&时域采样

我们在频域对信号进行采样之后，要能够进行信号的恢复，至于对于频域采样的信号能不能恢复成一个时域的信号，时域信号会不会发生混叠，这些就需要看频域采样定理。在学习频域的采样定理前，首先来回时域的采样定理——奈奎斯特定理。

### 2.1、时域采样定理

![3_时域采样定理](https://raw.githubusercontent.com/FunPlus007/Markdown4Zhihu/master/Data/DSP_3/3_时域采样定理.jpg)

所谓时域采样就是利用采样脉冲序p(t)从连续时间信号ft()中抽取一系列离散样本值的过程。这样得到的离散信号称为采样信号 <img src="https://www.zhihu.com/equation?tex=f_s(t)" alt="f_s(t)" class="ee_img tr_noresize" eeimg="1"> 。 <img src="https://www.zhihu.com/equation?tex=p_{Ts}" alt="p_{Ts}" class="ee_img tr_noresize" eeimg="1"> 为周期矩形采样信号（也可用周期冲激信号进行采样，周期冲激信号的傅里叶级数为1，为了**体现傅里叶级数的加权作用这里来看矩形信号**），下面三个式子为求周期信号的傅里叶变换。 <img src="https://www.zhihu.com/equation?tex=F_k" alt="F_k" class="ee_img tr_noresize" eeimg="1"> 为傅里叶系数。（周期信号的傅里叶变换由无穷多个冲激函数组成，这些冲激函数位于周期信号的各谐波频率 <img src="https://www.zhihu.com/equation?tex=nw_0" alt="nw_0" class="ee_img tr_noresize" eeimg="1"> (n=0,±1,±2,⋯)处，其强度对应信号展开为傅里叶级数）。


<img src="https://www.zhihu.com/equation?tex=\begin{align}
f_s(t) &= f(t)p_{rs}(t)\\
p_{Ts}(t) &= \sum_{n=-\infty}^{+\infty}g_\tau (t-nT_s)\\
P_{Ts}(jw) &= 2\pi \sum_{k=-\infty}^{+\infty}F_k\delta(w-kw_s)\\
F_k &= \frac{1}{T_s}\int_{-\frac{T_s}{2}}^{\frac{T_s}{2}} P_{Ts}(t)e^{jkw_st}dt\\
\end{align}
" alt="\begin{align}
f_s(t) &= f(t)p_{rs}(t)\\
p_{Ts}(t) &= \sum_{n=-\infty}^{+\infty}g_\tau (t-nT_s)\\
P_{Ts}(jw) &= 2\pi \sum_{k=-\infty}^{+\infty}F_k\delta(w-kw_s)\\
F_k &= \frac{1}{T_s}\int_{-\frac{T_s}{2}}^{\frac{T_s}{2}} P_{Ts}(t)e^{jkw_st}dt\\
\end{align}
" class="ee_img tr_noresize" eeimg="1">
根据上面的式子，时域的相乘为频域的卷积，可以从频域的角度分析这个采样信号 <img src="https://www.zhihu.com/equation?tex=F_s(jw)" alt="F_s(jw)" class="ee_img tr_noresize" eeimg="1"> 。（因为 <img src="https://www.zhihu.com/equation?tex=\delta" alt="\delta" class="ee_img tr_noresize" eeimg="1"> 函数的存在，卷积运算很简单）

<img src="https://www.zhihu.com/equation?tex=F_s(jw) &=& FT[f(t)p_{Ts}(t)] = \frac{1}{2\pi}F(jw)*P_{Ts}(jw)\\
&=&\frac{1}{2\pi} F(jw)*2\pi \sum_{k=-\infty}^{+\infty}F_k\delta(w-kw_s)\\
&=& \sum_{k=-\infty}^{+\infty}F_kF[j(w-kw_s)]
" alt="F_s(jw) &=& FT[f(t)p_{Ts}(t)] = \frac{1}{2\pi}F(jw)*P_{Ts}(jw)\\
&=&\frac{1}{2\pi} F(jw)*2\pi \sum_{k=-\infty}^{+\infty}F_k\delta(w-kw_s)\\
&=& \sum_{k=-\infty}^{+\infty}F_kF[j(w-kw_s)]
" class="ee_img tr_noresize" eeimg="1">
上式表明，连续信号 <img src="https://www.zhihu.com/equation?tex=f(t)" alt="f(t)" class="ee_img tr_noresize" eeimg="1"> 在被采样后，采样信号的频谱是连续信号的频 谱 <img src="https://www.zhihu.com/equation?tex=F(jw)" alt="F(jw)" class="ee_img tr_noresize" eeimg="1"> 以 <img src="https://www.zhihu.com/equation?tex=w_s" alt="w_s" class="ee_img tr_noresize" eeimg="1"> 为间隔进行**周期延拓**，在重复过程中被采样采样信号的傅里叶级数的系数所加权。

进行周期延拓，为了避免两个延拓的周期的混叠，就需要设计这个采样频率 <img src="https://www.zhihu.com/equation?tex=w_s" alt="w_s" class="ee_img tr_noresize" eeimg="1"> 。变得到上面的奈奎斯特采样频率。

下图是一个冲激信号对一个带限信号的时域采样示意图。

![3_冲击信号采样](https://raw.githubusercontent.com/FunPlus007/Markdown4Zhihu/master/Data/DSP_3/3_冲击信号采样.jpg)

时域采样信号的恢复过程如下图所示

![3_时域采样恢复](DSP_3/3_时域采样恢复-17206003686921.jpg)

可以发现关键在于设计滤波器（上面图示采用理想滤波器），滤除周期延拓的信号，然后进行逆变换。但是理想滤波器是物理不可实现的，所以就需要等到后面再介绍滤波器的设计。DSP的无非就两件事情，DFT/IDFT（FFT/IFFT）和数字滤波器设计（IIR、FIR）。

### 2.2、频域采样定理

简单回顾了一下信号的采样和恢复，现在让我们回到DFT这个东西来，前面也说到DFT的物理意义是对离散信号在单位圆（虚轴[0,2pi]）上的N点等间隔采样。

对于长度为M的离散序列x(n),当频域采样点数 <img src="https://www.zhihu.com/equation?tex=N>M" alt="N>M" class="ee_img tr_noresize" eeimg="1"> 时，即N点的DFT得到 <img src="https://www.zhihu.com/equation?tex=X(k)" alt="X(k)" class="ee_img tr_noresize" eeimg="1"> ，进行IDFT可以回复原序列x(n)，否则产生混叠。

<img src="https://www.zhihu.com/equation?tex=x_N(n)=IDFT[X(k)] = x(n)
" alt="x_N(n)=IDFT[X(k)] = x(n)
" class="ee_img tr_noresize" eeimg="1">
 推导参见下面的教材。






<img src="https://www.zhihu.com/equation?tex=X(k)$ <img src="https://www.zhihu.com/equation?tex=恢复成为" alt="恢复成为" class="ee_img tr_noresize" eeimg="1"> X(e^{jw})$，下面这个内插公式和FIR滤波器到的结构设计有很大关系。
" alt="X(k)$ <img src="https://www.zhihu.com/equation?tex=恢复成为" alt="恢复成为" class="ee_img tr_noresize" eeimg="1"> X(e^{jw})$，下面这个内插公式和FIR滤波器到的结构设计有很大关系。
" class="ee_img tr_noresize" eeimg="1">
\begin{align}
X(e^{jw}) = \sum_{k=0}^{N-1} X(k)\varphi(w-\frac{2\pi}{N}k)\\
\varphi(w) = \frac{1}{N} \frac{sin(wN/2)}{sin(w/2)} e^{-jw(\frac{N-1}{2})} \\

\end{align}

<img src="https://www.zhihu.com/equation?tex=## 3、DFT在信号谱分析中的应用



### 3.1 对连续模拟信号的谱分析

信号的谱分析就是计算信号的傅里叶变换。对于连续信号，进行时域采样进行离散之后，对其进行DFT，进行近似谱分析，使用DFT表示是出FT、DTFT、ZT。为什么要用DFT进行表示，一个是FFT的快速计算，另一个就是他是离散的。

下面注意 对离散的频域使用 <img src="https://www.zhihu.com/equation?tex=\omega" alt="\omega" class="ee_img tr_noresize" eeimg="1"> ,连续的频率使用 <img src="https://www.zhihu.com/equation?tex=\Omega" alt="\Omega" class="ee_img tr_noresize" eeimg="1"> 

为什么说是近似？

时域信号有限长，则频谱无限宽。时域信号无限长，频谱才能也有限宽。（比如一个cos进行FT之后，频谱就是两个冲激函数）。

对连续信号 <img src="https://www.zhihu.com/equation?tex=x_a(t)" alt="x_a(t)" class="ee_img tr_noresize" eeimg="1"> 进行采样，得到离散信号 <img src="https://www.zhihu.com/equation?tex=x(n) = x_a(nT)" alt="x(n) = x_a(nT)" class="ee_img tr_noresize" eeimg="1"> ，T是采样间隔。再对 <img src="https://www.zhihu.com/equation?tex=x(n)" alt="x(n)" class="ee_img tr_noresize" eeimg="1"> 进行DFT得到 <img src="https://www.zhihu.com/equation?tex=X(k)" alt="X(k)" class="ee_img tr_noresize" eeimg="1"> 。这里的 <img src="https://www.zhihu.com/equation?tex=X(k)" alt="X(k)" class="ee_img tr_noresize" eeimg="1"> 和 <img src="https://www.zhihu.com/equation?tex=x(n)" alt="x(n)" class="ee_img tr_noresize" eeimg="1"> 是有有限长的，但是根据上述理论分析， 信号时域和频域带宽的互补关系，要得到有限宽的频谱，时域信号要无限长，序列无限长之后，便无法满足DFT的变换条件。因此 DFT进行谱分析必然是近似的。要滤去原始模拟信号的的幅度很小的高频信号和幅度很小的部分时域信号。

由DFT恢复出FT，直接分析 <img src="https://www.zhihu.com/equation?tex=X(k)" alt="X(k)" class="ee_img tr_noresize" eeimg="1"> 是无法看到 <img src="https://www.zhihu.com/equation?tex=X_a(j \Omega)" alt="X_a(j \Omega)" class="ee_img tr_noresize" eeimg="1"> 的全部频谱特性的,需要乘上一个采样周期T，否则只能看到N个离散的采样点的谱线，这个是所谓的栅栏效应

同时对无限长的时域模拟信号进行上述的采样,DFT,恢复FT，是需要先对无限长信号进行截断，这会产生误差，这就是截断效应。

下面来分析这两个。

**栅栏效应**

由时域离散信号和模拟信号傅里叶变换之间的关系可得,T为采样时间间隔，模拟信号持续时间为 <img src="https://www.zhihu.com/equation?tex=T_p" alt="T_p" class="ee_img tr_noresize" eeimg="1"> ，那么采样点数为 <img src="https://www.zhihu.com/equation?tex=N = \frac{T_p}{T}" alt="N = \frac{T_p}{T}" class="ee_img tr_noresize" eeimg="1"> ：
" alt="## 3、DFT在信号谱分析中的应用



### 3.1 对连续模拟信号的谱分析

信号的谱分析就是计算信号的傅里叶变换。对于连续信号，进行时域采样进行离散之后，对其进行DFT，进行近似谱分析，使用DFT表示是出FT、DTFT、ZT。为什么要用DFT进行表示，一个是FFT的快速计算，另一个就是他是离散的。

下面注意 对离散的频域使用 <img src="https://www.zhihu.com/equation?tex=\omega" alt="\omega" class="ee_img tr_noresize" eeimg="1"> ,连续的频率使用 <img src="https://www.zhihu.com/equation?tex=\Omega" alt="\Omega" class="ee_img tr_noresize" eeimg="1"> 

为什么说是近似？

时域信号有限长，则频谱无限宽。时域信号无限长，频谱才能也有限宽。（比如一个cos进行FT之后，频谱就是两个冲激函数）。

对连续信号 <img src="https://www.zhihu.com/equation?tex=x_a(t)" alt="x_a(t)" class="ee_img tr_noresize" eeimg="1"> 进行采样，得到离散信号 <img src="https://www.zhihu.com/equation?tex=x(n) = x_a(nT)" alt="x(n) = x_a(nT)" class="ee_img tr_noresize" eeimg="1"> ，T是采样间隔。再对 <img src="https://www.zhihu.com/equation?tex=x(n)" alt="x(n)" class="ee_img tr_noresize" eeimg="1"> 进行DFT得到 <img src="https://www.zhihu.com/equation?tex=X(k)" alt="X(k)" class="ee_img tr_noresize" eeimg="1"> 。这里的 <img src="https://www.zhihu.com/equation?tex=X(k)" alt="X(k)" class="ee_img tr_noresize" eeimg="1"> 和 <img src="https://www.zhihu.com/equation?tex=x(n)" alt="x(n)" class="ee_img tr_noresize" eeimg="1"> 是有有限长的，但是根据上述理论分析， 信号时域和频域带宽的互补关系，要得到有限宽的频谱，时域信号要无限长，序列无限长之后，便无法满足DFT的变换条件。因此 DFT进行谱分析必然是近似的。要滤去原始模拟信号的的幅度很小的高频信号和幅度很小的部分时域信号。

由DFT恢复出FT，直接分析 <img src="https://www.zhihu.com/equation?tex=X(k)" alt="X(k)" class="ee_img tr_noresize" eeimg="1"> 是无法看到 <img src="https://www.zhihu.com/equation?tex=X_a(j \Omega)" alt="X_a(j \Omega)" class="ee_img tr_noresize" eeimg="1"> 的全部频谱特性的,需要乘上一个采样周期T，否则只能看到N个离散的采样点的谱线，这个是所谓的栅栏效应

同时对无限长的时域模拟信号进行上述的采样,DFT,恢复FT，是需要先对无限长信号进行截断，这会产生误差，这就是截断效应。

下面来分析这两个。

**栅栏效应**

由时域离散信号和模拟信号傅里叶变换之间的关系可得,T为采样时间间隔，模拟信号持续时间为 <img src="https://www.zhihu.com/equation?tex=T_p" alt="T_p" class="ee_img tr_noresize" eeimg="1"> ，那么采样点数为 <img src="https://www.zhihu.com/equation?tex=N = \frac{T_p}{T}" alt="N = \frac{T_p}{T}" class="ee_img tr_noresize" eeimg="1"> ：
" class="ee_img tr_noresize" eeimg="1">
X(e^{jw}) = \frac{1}{T}\sum_{m=-\infty}^{+\infty}X_a[j(\frac{w}{T} - \frac{2\pi m}{T})]

<img src="https://www.zhihu.com/equation?tex= <img src="https://www.zhihu.com/equation?tex=w = \Omega T" alt="w = \Omega T" class="ee_img tr_noresize" eeimg="1"> 
" alt=" <img src="https://www.zhihu.com/equation?tex=w = \Omega T" alt="w = \Omega T" class="ee_img tr_noresize" eeimg="1"> 
" class="ee_img tr_noresize" eeimg="1">
\begin{align}
X(k) =& DFT[X(n)]_N = X(e^{jw})|_{w=2\pi k/N} = X(e^{j\frac{2\pi k}{N}})\\
=& \frac{1}{T}\sum_{m=-\infty}^{+\infty}X_a[j(\frac{2\pi k/N}{T} - \frac{2\pi m}{T})\\
=& \frac{1}{T}\hat{X}_a(\frac{2\pi}{T_p}k)

\end{align}
$$
所以可以看出，对模型信号进行采样后在进行DFT，要恢复原始模拟信号的FT，需要将X(k)乘上采样间隔T，可以得到模拟信号频谱周期延拓的第一个周期，如下图所示，否则无法得到原始模拟信号的的全部频谱，只有N个离散的谱线。



**截![dsp_3_DFT分析连续信号频谱](https://raw.githubusercontent.com/FunPlus007/Markdown4Zhihu/master/Data/DSP_3/dsp_3_DFT分析连续信号频谱.jpg)断效应**

对一个无限长的信号继续截断之后采样DFT，必然会产生误差，这就是截断效应，避免截断效应一个是加长信号的持续时间，让截断的信号尽可能和原始信号相近，第二种便是加窗，截断是用一个矩形窗去截，在边界，矩形窗是之间垂直下降的，这样变化太突兀了，可以设计窗函数，就是设计截断窗的形状，优化在边界的滚降特性，减小误差。

![3_截断效应](https://raw.githubusercontent.com/FunPlus007/Markdown4Zhihu/master/Data/DSP_3/3_截断效应.jpg)

一个Sa函数的FT理论上来说是一个门函数，但是Sa是无限长的，计算机算不了的，所以将Sa截断，用的矩形窗截断（乘的一个R），如上面的（a）所示,就是截断后的Sa函数，然后对它进行采样，DFT，乘上T（采样间隔）恢复模拟信号的频谱。可以看到（c）就是计算机应该操作出来的结果，（b）是理论上的频谱。

可以发现，他们两个的区别还是很大的，高频整个部分都有波动。改善的办法就是设计窗函数，经典的升余弦窗（汉宁Hanning窗），改进的升余弦窗（Hamming窗），布莱克曼窗，贝塞尔窗等等。这些详细介绍在窗函数法设计FIR滤波器。

至此，我们以及分析出了如何让计算机处理一个模拟信号，对他进行采样之后，进行DFT，进而可以恢复他的FT。选择合适的采样频率，避免栅栏效应和截断效应减小误差。

