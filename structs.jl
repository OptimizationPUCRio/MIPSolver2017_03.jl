mutable struct prob_lp
    A
    c
    n
    m
    xlb
    xub
    l
    u
end

mutable struct resposta_relaxado
    obj
    vars
    status
end

mutable struct modelo_lista
    problem::prob_lp
    resp::resposta_relaxado
    modelo_lista() = new()
end
