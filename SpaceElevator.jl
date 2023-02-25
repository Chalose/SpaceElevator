#= SpaceElveator的粗略计算
假设：
1、电梯为抗拉型巨构；
2、货箱离地高度忽略为0；
3、线缆均相材质，不同轨道位置r其截面积S(r), 有S(r)=T(r)/P;    T(r)为线缆张力，P为材料抗拉极限
4、BC段配重完全由线缆质量提供，使得电梯质心维持在同步轨道Rₕ处
    C 。最远点
      |
      |
      |
      |
      |
      |
    B 。同步轨道锚点
      |
      |
      |
      |
      |
    A 。近地点
目的：
i、main1()对比了AB段钢材与CNT的性能差异(线径变化、轴向张力)
ii、main2()分析了全局CNT线缆的轴向张力，以及采用不同强度CNT的消耗的质量差异
=#
ENV["GKS_ENCODING"] = "utf-8"
using Plots
using Roots
plotlyjs()
include("IntGL.jl")

# 全局常数(SI)
R₀ = 6.7310e6    # 地球平均半径
Rₕ = 4.21640e7    # 地球同步轨道半径
G = 6.67259e-11  # 引力常数
M = 5.965e24     # 地球质量
P₁ = 600e6       # 钢材抗拉强度
P₂ = 25*P₁       # 碳纳米管抗拉强度
ρ₁ = 7.850e3     # 钢材平均密度
ρ₂ = 1/6*ρ₁      # 碳纳米管密度
ω = 7.292e-5     # 自转角速度
MA =20e3;  F₀ = MA*9.8  # 货箱质量以及对线缆负载
MB = 1000e3       # 同步轨道锚点质量

α₁ = G*M*ρ₁/P₁   # 常数α,β: 下标1为钢，下标2为CNT
α₂ = G*M*ρ₂/P₂
β₁ = ω^2*ρ₁/P₁
β₂ = ω^2*ρ₂/P₂

# AB段：CNT与钢材性能对比
function main1()
    # 求解截面积函数 S(r)、最大张力函数 T(r)
    dr = 100e3
    h = (Vector{Float64}(R₀:dr:Rₕ) .- R₀)./1000   # 离地高度h (km)
    
    S₁ = [F₀/P₁ * exp(α₁*(1/R₀ - 1/r) + β₁/2*(R₀^2 - r^2)) for r = R₀:dr:Rₕ]
    S₂ = [F₀/P₂ * exp(α₂*(1/R₀ - 1/r) + β₂/2*(R₀^2 - r^2)) for r = R₀:dr:Rₕ]
    T₁ = P₁*S₁
    T₂ = P₂*S₂
    
    f1(r) = F₀/P₁ * exp.(α₁*(-1/r .+ 1/R₀) .+ β₁/2*(-r.^2 .+ R₀^2)) * ρ₁
    f2(r) = F₀/P₂ * exp.(α₂*(-1/r .+ 1/R₀) .+ β₂/2*(-r.^2 .+ R₀^2)) * ρ₂

    # 钢与碳纳米管对比图
    p1 = plot(h, 2*sqrt.(S₁/π), 
        title="线径", 
        ylabel="Φ (m)",  
        label="Steel")
    p2 = plot(h, 2*sqrt.(S₂/π)*100,
        xlabel="h (km)", 
        ylabel="Φ (cm)", 
        color="red",  
        label="CNT")
    p3 = plot(h, T₁, 
        title="轴向张力", 
        ylabel="T (N)", 
        label="Steel")
    p4 = plot(h, T₂, 
        xlabel="h (km)", 
        ylabel="T (N)", 
        color="red", 
        label="CNT")
    fig1 = plot(p1, p3, p2, p4, layout=(2, 2), legend=:outertopright)
    savefig(fig1, "SpaceElevator01.png")
end


# 考虑平衡配重后的CNT电梯
function main2()
    dr = 100e3

    T_AB(r) = F₀*exp(α₂*(1/R₀ - 1/r) + β₂/2*(R₀^2 - r^2))
    F₁ = MB*G*M/Rₕ^2
    F₂ = T_AB(Rₕ)
    T_BC(r) = (F₁+F₂)*exp(α₂*(1/Rₕ - 1/r) + β₂/2*(Rₕ^2 - r^2))

    # 确定最远点Rc位置
    S_AB(r) = T_AB(r)/P₂
    S_BC(r) = T_BC(r)/P₂
    fun1(r) = (r - Rₕ)*S_AB(r)
    c1 = IntGL(fun1, R₀, Rₕ, 400)
    fun2(r) = (r - Rₕ)*S_BC(r)
    f_Rc(x) = c1 + IntGL(fun2, Rₕ, x, 400)
    Rc = find_zero(f_Rc, 1.5*Rₕ)
    println("电梯最远端离地 ", Rc/1e3, " (km)")

    # 计算全局张力(由于MB的作用，张力曲线有不连续性)
    T1 = [T_AB(r) for r = R₀:dr:Rₕ]
    T2 = [T_BC(r) for r = Rₕ+dr:dr:Rc]
    T = [T1; T2]
    lt = length(T)    # 避免数组维度问题
    h = [(0.0 + (i-1)*dr)/1e3 for i=1:lt]
    p1 = plot(h, T,
        xlabel="h (km)", 
        ylabel="T (N)",
        color="red", 
        label="25倍于钢强度时的张力曲线")
    
    # 计算不同CNT强度时,消耗线缆总质量m
    k = 25:50     # 倍率
    i = 1
    m = zeros(length(k), 1)
    for j in k
        P₂ = j*P₁
        α₂ = G*M*ρ₂/P₂
        β₂ = ω^2*ρ₂/P₂

        fk1(r) = F₀*exp(α₂*(1/R₀ - 1/r) + β₂/2*(R₀^2 - r^2)) / P₂ * ρ₂
        fk2(r) = (F₁+F₂)*exp(α₂*(1/Rₕ - 1/r) + β₂/2*(Rₕ^2 - r^2)) / P₂ * ρ₂

        m[i] = ((IntGL(fk1, R₀, Rₕ, 400))  + IntGL(fk2, Rₕ, Rc, 400))/1e3
        i += 1
    end
    p2 = plot(k, m,
        xlabel="k", 
        ylabel="m (t)",
        color="red",
        label="k倍于钢强度时的线缆质量")

    fig2 = plot(p1, p2, 
        layout=(2, 1), 
        legend=:outertopright)
    savefig(fig2, "SpaceElevator02.png")
end
main1()
main2()
