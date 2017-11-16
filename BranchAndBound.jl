# Branch and Bound, Eduardo Brito.
module SolverBrito

    using JuMP, Gurobi
#--------------------------------------------------------------------------------------------------
    # Struct para informacoes de um modelo generico
    mutable struct prob_lp
        A
        c
        n
        m
        xlb
        xub
        rowlb
        rowub
        vtypes
        solver
    end

#--------------------------------------------------------------------------------------------------
    # Struct para resposta de um problema generico
    mutable struct resposta_relaxado
        obj #Z*
        vars #X*
        status #Status, :Optimal, :Infeasible, ...
    end

#--------------------------------------------------------------------------------------------------
    # Struct que combina o modelo com a resposta
    mutable struct modelo_lista
        problem::prob_lp
        resp::resposta_relaxado
        modelo_lista() = new()
    end

#--------------------------------------------------------------------------------------------------
    # Funcao para extrair informacoes do modelo e passar para a struct prob_lp
    function converte_modelo(m::Model)
        md = deepcopy(m)

        if md.objSense == :Max
            @objective(md,:Min,-md.obj)
        end

        c = JuMP.prepAffObjective(md)
        A = JuMP.prepConstrMatrix(md)
        n, m = size(A)
        xlb = copy(md.colLower)
        xub = copy(md.colUpper)
        rowlb, rowub = JuMP.prepConstrBounds(md)
        solver = md.solver
        vtypes = md.colCat .!= :Cont
        problem = prob_lp(A,c,n,m,xlb,xub,rowlb,rowub,vtypes,solver)
        return problem
    end

#--------------------------------------------------------------------------------------------------
    # Funcao que resolve um modelo do struct "prob_lp", e devolve em forma do struct resp
    function solve_relax(modelo::prob_lp)
        model = Model(solver=modelo.solver)
        @variable(model, x[1:modelo.m])
        for i in 1:modelo.m
            setlowerbound(x[i], modelo.xlb[i])
            setupperbound(x[i], modelo.xub[i])
        end
        @constraint(model, modelo.rowlb .<= modelo.A*x .<= modelo.rowub)
        @objective(model, Min, dot(modelo.c,x))

        status = solve(model)
        resposta = resposta_relaxado(getobjectivevalue(model),getvalue(x),status)
        return resposta
    end

#--------------------------------------------------------------------------------------------------
    # Funcao podas. Recebe a resposta do filho e pode por inviabilidade, limite e otimalidade.
    function podas(resp ::resposta_relaxado, global_bound,vtypes)

        #Poda por Inviabilidade
        if(resp.status != :Optimal)
            return "erro"
        end
        #Poda por Limite
        if (resp.obj>global_bound[2])
            return "erro"
        end
        #Poda por Otimalidade
        if ( sum(abs.(resp.vars[vtypes] - round.(resp.vars[vtypes]))) == 0)
            return Float64(resp.obj)
        end
        return "sucesso"
    end

#--------------------------------------------------------------------------------------------------
    # Funcao exporta_model. Recebe o modelo original e atualiza os valores do z*, x, ... , e preenche listas dos nos e solucoes inteiras
    function exporta_model(m,zinf,nodes,integer_solutions,time,global_bound,iter)
        m.objVal = zinf.resp.obj
        m.colVal = zinf.resp.vars
        m.objBound = abs(global_bound[1]


        m.ext[:status] = zinf.resp.status
        m.ext[:time] = time
        m.ext[:nodes] = size(nodes)[1]
        m.ext[:intsols] = size(integer_solutions)[1]
        m.ext[:iter] = iter

        return m
    end
#--------------------------------------------------------------------------------------------------
    # Funcao que recebe a lista de problemas e retorna o indice do problema a ser trabalhado. (Para P2)
    function acha_problema(lista)
        return 1
    end

#--------------------------------------------------------------------------------------------------
    # Funcao que recebe a lista de problemas e o indice do problema selecionado e retorna o indice da variavel para o branch. (Pega a mais fracionaria)
    function acha_variavel(lista,ind_prob)
        X = lista[ind_prob].resp.vars
        Var_Cont = convert(Array{Int,1},lista[ind_prob].problem.vtypes)

        ind = indmax(abs.(X.*Var_Cont - round.(X.*Var_Cont)))
        return ind
    end

#--------------------------------------------------------------------------------------------------
    #Funcao SOLVEMIP
    function SolveMIP(m)
        tic()

        # Cria listas: De problemas, nos encontrados e solucoes inteiras. -------------------------
        lista = Vector{modelo_lista}()
        nodes = Vector{modelo_lista}()
        integer_solutions = Vector{modelo_lista}()

        # Adiciona o 1o problema de forma manual na lista. ----------------------------------------
        zinf = modelo_lista()
        zinf.problem = converte_modelo(m)
        zinf.resp = solve_relax(zinf.problem)
        global_bound = [zinf.resp.obj,+Inf]
        push!(lista,zinf)

        # Comeca o algoritmo do branch and bound. -------------------------------------------------
        iter = 0
        if sum(zinf.problem.vtypes) != 0 #Esta checando se todas as variaveis ja sao inteiras

            while (abs(global_bound[2] - global_bound[1]) >= exp10(-5) && size(lista)[1] != 0 && iter <= 1000) #Para por limite de iteracoes, bound e tamanho da lista

                # Seleciona problema usando a funcao "acha_problema". -----------------------------
                ind_prob = acha_problema(lista)

                # Seleciona a variavel a ser modificada -------------------------------------------
                ind = acha_variavel(lista,ind_prob)

                # Branch --------------------------------------------------------------------------
                prob_LF = deepcopy(lista[ind_prob].problem)
                prob_LF.xlb[ind] = floor(lista[ind_prob].resp.vars[ind]) + 1 # Altera bounds da variavel

                prob_RT = deepcopy(lista[ind_prob].problem)
                prob_RT.xub[ind] = floor(lista[ind_prob].resp.vars[ind]) # Altera bounds da variavel

                # Soluciona os filhos ------------------------------------------------------------
                resp_LF = solve_relax(prob_LF)
                resp_RT = solve_relax(prob_RT)

                # Podas --------------------------------------------------------------------------
                poda_LF = podas(resp_LF, global_bound,prob_LF.vtypes)
                poda_RT = podas(resp_RT, global_bound,prob_RT.vtypes)

                # Monta novas variaveis da struct modelo_lista, prontas para serem adicionadas na lista
                #LEFT
                novo_LF = modelo_lista()
                novo_LF.problem = prob_LF
                novo_LF.resp = resp_LF

                #RIGHT
                novo_RT = modelo_lista()
                novo_RT.problem = prob_RT
                novo_RT.resp = resp_RT

                # Olha para as podas e atualiza o global_bound -----------------------------------
                if typeof(poda_RT) == Float64
                    if typeof(poda_LF) == Float64
                        if poda_RT <= poda_LF
                            global_bound[2] = poda_RT
                        else
                            global_bound[2] = poda_LF
                        end
                    else
                        global_bound[2] = poda_RT
                    end
                elseif typeof(poda_LF) == Float64
                        global_bound[2] = poda_LF
                end

                #Atualiza Zinf (Melhor solucao ate o momento) ----------------------------------
                if typeof(poda_LF) == Float64
                    if novo_LF.resp.obj <= global_bound[2]
                        zinf = deepcopy(novo_LF)
                    end
                end
                if typeof(poda_RT) == Float64
                    if typeof(poda_LF) == Float64
                        if novo_RT.resp.obj <= global_bound[2] && novo_RT.resp.obj <= novo_LF.resp.obj
                            zinf = deepcopy(novo_RT)
                        end
                    else
                        if novo_RT.resp.obj <= global_bound[2]
                            zinf = deepcopy(novo_RT)
                        end
                    end
                end

                #Remove o problema original e adiciona os novos ---------------------------------
                # remove
                push!(nodes,lista[ind_prob])
                deleteat!(lista,ind_prob)

                # adiciona
                if poda_LF == "sucesso"
                    push!(lista,novo_LF)
                elseif typeof(poda_LF) == Float64
                    push!(integer_solutions, novo_LF)
                end

                if poda_RT == "sucesso"
                    push!(lista,novo_RT)
                elseif typeof(poda_RT) == Float64
                    push!(integer_solutions, novo_RT)
                end
                iter += 1
            end
        end

        # Fim do while, agora verificar model.sense, corrigir caso necessario o z*.
        if m.objSense == :Max
            zinf.resp.obj = - zinf.resp.obj
        end

        time = toc() # Registro do tempo

        # Usa a funcao exporta_model para compor o modelo final usando a melhor solucao encontrada ate agora.
        model = exporta_model(m,zinf,nodes,integer_solutions,time,global_bound,iter)

        return model
    end
end
