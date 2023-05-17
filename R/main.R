
flux_corr_main <- function(flux_solutions) {
  ## func1
  # n = as.data.frame(t(flux_solutions))
  # tmp <- cor(n)
  # tmp[upper.tri(tmp)] <- 0
  # diag(tmp) <- 0
  # res_corr <- n[,!apply(tmp,2,function(x) any(abs(x) > 0.3))]

  ## func2
  library('caret')
  n = as.data.frame(t(flux_solutions))
  df2 <- cor(n)
  # hc = findCorrelation(df2, cutoff=0.3) # putt any value as a "cutoff"
  # hc = findCorrelation(df2, cutoff=0.3, exact = TRUE)

  hc = findCorrelation(df2, cutoff = 0.3) # putt any value as a "cutoff"
  hc = findCorrelation(df2, cutoff = 0.3, exact = TRUE) # putt any value as a "cutoff"
  hc = sort(hc)
  s_names = colnames(n)
  s_names = s_names[-c(hc)]
  return(s_names)
}

build_net <- function(net, anchor_rid) {
  all_compounds = rownames(net)
  #all_compounds = c(all_compounds, "upstream", "downstream")
  all_compounds = c(all_compounds, "downstream")
  ## vertex
  library(igraph)
  g <- make_empty_graph(n = length(all_compounds))
  g <- set.vertex.attribute(g, "name", value = all_compounds)
  V(g)$color <- "yellow"

  all_reactions = colnames(net)
  cm_inout = data.frame()
  ## edges
  for (i in 1:length(all_reactions)) {
    rid = all_reactions[i]

    c_in_nodes = all_compounds[which(net[, i] == -1)]
    c_out_nodes = all_compounds[which(net[, i] == 1)]
    if (identical(c_in_nodes, character(0))) {
      #upstream
      c_in_nodes = "upstream"
    }
    if (identical(c_out_nodes, character(0))) {
      #downstream
      c_out_nodes = "downstream"
    }
    for (c_in in c_in_nodes) {
      for (c_out in c_out_nodes) {
        #print(paste0(c_in, "->", c_out))
        cm_row = data.frame()
        cm_row[1, "C_in"] = c_in
        cm_row[1, "C_out"] = c_out
        cm_inout = rbind(cm_inout, cm_row)

        if (!are.connected(g, c_in, c_out)) {
          g <- add_edges(g, c(c_in, c_out), rid = rid)
        } else {
          old_rid = unlist(strsplit(E(g)[get.edge.ids(g, c(c_in, c_out))]$rid, split = ";"))
          new_rid = union(old_rid, rid)
          E(g)[get.edge.ids(g, c(c_in, c_out))]$rid = paste(new_rid, collapse = ";")
        }
      }
    }
  }

  ## anchor
  E(g)$color = "grey"
  E(g)[rid %in% anchor_rid]$color = "Salmon"

  up_nodes = neighbors(g, V(g)[name == "upstream"], mode = "out")$name
  down_nodes = neighbors(g, V(g)[name == "downstream"], mode = "in")$name
  up_down_nodes = intersect(up_nodes, down_nodes)

  V(g)[name %in% up_nodes]$color <- "white"
  V(g)[name %in% down_nodes]$color <- "Coral"
  V(g)[name %in% up_down_nodes]$color <- "Wheat"

  V(g)[name == "upstream"]$color <- "pink"
  V(g)[name == "downstream"]$color <- "grey"

  for (i in 1:length(rownames(cm_inout))) {
    cm_inout[i, ""]
  }

  #cm_inout = cm_inout[-which((cm_inout$C_in %in% c("upstream","downstream")) | (cm_inout$C_out %in% c("upstream","downstream"))),]

  res = list()
  res[["g"]] = g
  res[["cm"]] = cm_inout

  return(res)
}




find_anchor_rid <- function(net) {
  ## anchor edge
  nozero_count = apply(net, 2, function(c)
    sum(c != 0))
  anchor_rid = names(nozero_count[which(nozero_count > 2)])
  return(anchor_rid)
}

find_upstream_rid <- function(net) {
  ## anchor edge
  nozero_count = apply(net, 2, function(c)
    sum(c != 0))
  srid = names(nozero_count[which(nozero_count == 1)])
  upstream_rid = srid[net[, srid][net[, srid] != 0] == 1]

  return(upstream_rid)
}

find_downstream_rid <- function(net) {
  ## anchor edge
  nozero_count = apply(net, 2, function(c)
    sum(c != 0))
  srid = names(nozero_count[which(nozero_count == 1)])
  downstream_rid = srid[net[, srid][net[, srid] != 0] == -1]

  return(downstream_rid)
}

add_upstream_row <- function(cmMat, upstream_rid) {
  cpds = rownames(cmMat)
  cmMat = rbind(cmMat, 0)
  cmMat[nrow(cmMat), upstream_rid] = -1
  row.names(cmMat) = c(cpds, "upstream")
  return(cmMat)
}


get_template_cm <- function(hsa_net, connected_compounds) {
  template_cm = matrix(
    0,
    nrow = length(connected_compounds),
    ncol = length(connected_compounds),
    dimnames = list(connected_compounds, connected_compounds)
  )

  connected_net = hsa_net[which((hsa_net$C_in %in% connected_compounds) &
                                  (hsa_net$C_out %in% connected_compounds)), ]
  for (i in 1:length(rownames(connected_net))) {
    temp_cin = as.character(connected_net[i, "C_in"])
    temp_cout = as.character(connected_net[i, "C_out"])

    template_cm[temp_cin, temp_cout] = 1
    template_cm[temp_cout, temp_cin] = 1
  }

  return(template_cm)

}

## cycle
get_anchor_cycle_affiliated_edges <-
  function(hsa_net, g, anchor_rid, anchor_node) {
    #anchor_node = unique(as.vector(ends(g, E(g)[rid%in%anchor_rid], names = TRUE)))

    library(igraph)
    compound_graph <- as.undirected(g)
    compound_components = igraph::groups(components(compound_graph))
    connected_compounds = as.vector(compound_components[[1]])

    template_cm = get_template_cm(hsa_net, connected_compounds)

    library('ggm')
    fund_cycles = fundCycles(template_cm)
    cycles = fund_cycles[lapply(fund_cycles, length) > 2]

    replace_cname <- function(x) {
      x <- connected_compounds[x]
      y = c(x[1:(length(x) / 2)], x[1])
      return(y)
    }
    cycles_cname <- lapply(cycles, replace_cname)

    # judge directed cycle
    judge_directed_cycle <- function(x) {
      pos_judge = TRUE
      neg_judge = TRUE

      for (i in 1:(length(x) - 1)) {
        edge_c = c(x[i], x[i + 1])
        ud_edge = hsa_net[which(((hsa_net$C_in %in% edge_c) &
                                   (hsa_net$C_out %in% edge_c))), ]

        d_edge = hsa_net[which(((hsa_net$C_in == x[i]) &
                                  (hsa_net$C_out == x[i + 1]))), ]
        if (!(TRUE %in% ud_edge$Irreversible) &
            length(rownames(d_edge)) == 0) {
          pos_judge = FALSE
        }
      }

      neg_x = rev(x)
      for (i in 1:(length(neg_x) - 1)) {
        edge_c = c(neg_x[i], neg_x[i + 1])
        ud_edge = hsa_net[which(((hsa_net$C_in %in% edge_c) &
                                   (hsa_net$C_out %in% edge_c))), ]

        d_edge = hsa_net[which(((hsa_net$C_in == neg_x[i]) &
                                  (hsa_net$C_out == neg_x[i + 1]))), ]
        if (!(TRUE %in% ud_edge$Irreversible) &
            length(rownames(d_edge)) == 0) {
          neg_judge = FALSE
        }

      }

      if (pos_judge | neg_judge) {
        return(x)
      }

    }
    cycles_cname_directed = lapply(cycles_cname, judge_directed_cycle)
    cycles_cname_directed = cycles_cname_directed[lapply(cycles_cname_directed, length) > 0]

    anchor_cycle_affiliated_edges = c()
    if (length(cycles_cname_directed) > 0) {
      for (i in 1:length(cycles_cname_directed)) {
        g1 = subgraph(g, V(g)[name %in% cycles_cname_directed[[i]]])

        if (length(intersect(E(g1)$rid, anchor_rid)) > 0) {
          for (cn in cycles_cname_directed[[i]]) {
            in_nodes = neighbors(g, V(g)[name == cn], mode = "in")$name
            in_nodes = setdiff(in_nodes,  cycles_cname_directed[[i]])
            in_nodes =  setdiff(in_nodes, anchor_node)

            for (ind in in_nodes) {
              in_suc = bfs(
                g,
                root = ind,
                "in",
                unreachable = FALSE,
                order = TRUE,
                rank = TRUE,
                father = TRUE,
                pred = TRUE,
                succ = TRUE,
                dist = TRUE
              )
              in_suc = na.omit(in_suc$order)
              in_suc = in_suc$name

              in_suc = c(cn, in_suc)
              g_in = subgraph(g, V(g)[name %in% in_suc])
              #plot(g_in)
              anchor_cycle_affiliated_edges =  union(E(g_in)$rid, anchor_cycle_affiliated_edges)
            }

            out_nodes = neighbors(g, V(g)[name == cn], mode = "out")$name
            out_nodes = setdiff(out_nodes,  cycles_cname_directed[[i]])
            out_nodes =  setdiff(out_nodes, anchor_node)

            for (od in out_nodes) {
              out_suc = bfs(
                g,
                root = out_nodes[1],
                "out",
                unreachable = FALSE,
                order = TRUE,
                rank = TRUE,
                father = TRUE,
                pred = TRUE,
                succ = TRUE,
                dist = TRUE
              )
              out_suc = na.omit(out_suc$order)
              out_suc = out_suc$name

              out_suc = c(cn, out_suc)
              g_out = subgraph(g, V(g)[name %in% out_suc])
              #plot(g_out)
              anchor_cycle_affiliated_edges =  union(E(g_out)$rid, anchor_cycle_affiliated_edges)
            }
          }
        }
      }
    }


    aca_edges = c()
    if (length(anchor_cycle_affiliated_edges) != 0) {
      for (i in 1:length(anchor_cycle_affiliated_edges)) {
        tmp_aca_edge = unlist(strsplit(anchor_cycle_affiliated_edges[i], split = ";"))
        aca_edges = union(aca_edges, tmp_aca_edge)
      }
    }

    return(aca_edges)


  }










# get_bvls_solution <- function(upstream_value, boundary_value, cmMat, anchor_cycle_affiliated_edges, error=FALSE) {
#   A = as.matrix(cmMat)
#   #upstream_value = 10
#
#   p=ncol(cmMat)
#   if(error) {
#     B = matrix(c(runif(nrow(cmMat)-1, min=0, max=0.000005),-upstream_value), nr=nrow(cmMat), nc=1)
#   }else {
#     B = matrix(c(rep(0,nrow(cmMat)-1),-upstream_value), nr=nrow(cmMat), nc=1)
#   }
#   boundary_low = rep(0,p)
#   free_edge_no = which(!colnames(A)%in%anchor_cycle_affiliated_edges)
#   boundary_low[free_edge_no] = boundary_value
#   #boundary_up = rep(Inf, p)
#   boundary_up = rep(upstream_value, p)
#
#   library(bvls)
#   bvls_varible <- bvls(A, B, bl = boundary_low, bu = boundary_up)
#
#   return(bvls_varible$x)
# }

get_bvls_solution <-
  function(upstream_value,
           boundary_value,
           cmMat,
           anchor_cycle_affiliated_edges,
           error = FALSE) {
    A = as.matrix(cmMat)
    #upstream_value = 10
    p = ncol(cmMat)
    # if(error) {
    #   B = matrix(c(runif(nrow(cmMat)-1, min=0, max=0.000005),-upstream_value), nr=nrow(cmMat), nc=1)
    # }
    B = matrix(c(rep(0, nrow(cmMat) - 1), -upstream_value),
               nr = nrow(cmMat),
               nc = 1)

    boundary_low = rep(0, p)
    free_edge_no = which(!colnames(A) %in% anchor_cycle_affiliated_edges)

    ## new tool
    boundary_low[free_edge_no] = boundary_value
    boundary_up = boundary_low + upstream_value

    library(bvls)
    bvls_varible <- bvls(A, B, bl = boundary_low, bu = boundary_up)

    return(bvls_varible$x)
  }



plot_flux_sim <-
  function(solutions, solutions_error, g, graph_name) {
    # g = delete_vertices(g, V(g)[name=="upstream"])
    # g = delete_vertices(g, V(g)[name=="downstream"])
    ge = g
    #graph_name = "flux_sim.png"
    for (i in 1:length(E(g)$rid)) {
      edge_name = E(g)$rid[i]
      #print(i)
      #print(edge_name)
      rids = unlist(strsplit(E(g)$rid[i], split = ";"))
      for (j in 1:length(rids)) {
        rid_no = as.numeric(unlist(strsplit(rids[j], split = "_"))[2])
        E(g)[rid == edge_name]$value = round(solutions[rid_no], 2)
        E(ge)[rid == edge_name]$value = round(solutions_error[rid_no], 2)
      }
    }

    library(igraph)
    library(qgraph)
    plot_name = file.path(paste0(graph_name, ".png"))
    png_size = ceiling(vcount(g) / 8)
    png_size[png_size < 15] = 15
    png(
      plot_name,
      height = png_size,
      width = png_size,
      units = "in",
      res = 200
    )
    e <- get.edgelist(g, names = FALSE)
    #l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),area=(18*vcount(g)^2),repulse.rad=(vcount(g)^3.2))
    #l = layout_nicely(g)
    #l = layout_with_kk(g)
    l = layout_with_graphopt(g)

    #plot(g, layout=l, vertex.size=2, edge.arrow.size=0.8, edge.width=E(g)$size)
    par(mfrow = c(1, 2))
    plot(
      g,
      edge.label = paste0(E(g)$rid, ":", E(g)$value),
      vertex.label = V(g)$name,
      vertex.size = 5,
      edge.arrow.size = 1,
      edge.width = 2,
      layout = l
    )
    mtext("origin", side = 1)
    plot(
      ge,
      edge.label = paste0(E(ge)$rid, ":", E(ge)$value),
      vertex.label = V(ge)$name,
      vertex.size = 5,
      edge.arrow.size = 1,
      edge.width = 2,
      layout = l
    )
    mtext("error", side = 1)

    dev.off()
  }


## write
write_flux_sim <- function(solutions, solutions_error, csv_name) {
  X = data.frame("flux" = solutions,
                 "flux_error" = solutions_error)
  rownames(X) = paste0("M_", c(1:ncol(cmMat)))
  write.table(
    X,
    paste0(csv_name, "_sim.csv"),
    row.names = TRUE,
    col.names = FALSE,
    sep = ","
  )
}

## main
library(readr)
flux_sim_main <-
  function(accurate_1,
           solution_num,
           metabolism_network_path) {
    #filename = unlist(strsplit(csvname, split = '.', fixed = T))[1] #'GSL3_cmMat'
    #cmMat <- read_csv(paste0(csvpath, csvname), col_names = TRUE)
    cmMat <- read_csv(metabolism_network_path, col_names = TRUE)
    get_cmMat <- function(GS_cmMat) {
      cmMat = as.data.frame(GS_cmMat)
      rownames(cmMat) = cmMat[, 1]
      cmMat = cmMat[, -1]
      return(cmMat)
    }
    cmMat = get_cmMat(cmMat)

    anchor_rid = find_anchor_rid(cmMat)
    upstream_rid = find_upstream_rid(cmMat)
    downstream_rid = find_downstream_rid(cmMat)

    cmMat = add_upstream_row(cmMat, upstream_rid)
    cmMat_res = build_net(cmMat, anchor_rid)
    g = cmMat_res[["g"]]
    cmMat_inout = cmMat_res[["cm"]]

    anchor_node = unique(as.vector(ends(g, E(g)[rid %in% anchor_rid], names = TRUE)))
    anchor_cycle_affiliated_edges = get_anchor_cycle_affiliated_edges(cmMat_inout, g, anchor_rid, anchor_node)

    times = solution_num
    anchor_num = length(anchor_rid)
    aca_edge_num = length(anchor_cycle_affiliated_edges)
    free_edge_num = ncol(cmMat) - aca_edge_num

    flux_df = data.frame(rep(0, ncol(cmMat)))
    flux_df = data.frame(colnames(cmMat))

    for (i in 1:times) {
      upstream_value = runif(1,
                             min = min(0.05 * free_edge_num, 1.0),
                             max = 1.5)
      #boundary_value = runif(free_edge_num, min=0.3, max=0.8)
      boundary_value = runif(free_edge_num, min = 0.00001, max = 0.000012)

      bvls_solution = get_bvls_solution(upstream_value,
                                        boundary_value,
                                        cmMat,
                                        anchor_cycle_affiliated_edges,
                                        FALSE)

      A = as.matrix(cmMat)
      B = matrix(c(rep(0, nrow(cmMat) - 1), -upstream_value),
                 nr = nrow(cmMat),
                 nc = 1)
      solution_bias = abs((A %*% bvls_solution) - B)
      if ((max(solution_bias) < accurate_1)) {
        flux_df = cbind(flux_df, bvls_solution)
      }
    }
    print(ncol(flux_df))

    rownames(flux_df) = flux_df[, 1]
    flux_df = flux_df[, -1]

    flux_df = as.data.frame(t(flux_df))
    rownames(flux_df) <- paste("S_", 1:nrow(flux_df), sep = "")

    wcorr_matrix_rowname = flux_corr_main(flux_df)
    wcorr_flux_df = flux_df[wcorr_matrix_rowname,]

    #write.table(wcorr_flux_df, paste0("simresult_flux/", filename, "_sim.csv"),row.names=FALSE,col.names=TRUE,sep=",")

    return(wcorr_flux_df)
  }

get_mod_genes_table <-
  function(module_gene_json,
           mod_genes_name,
           genes_grad) {
    mods_name = names(module_gene_json)
    mod_genes_table = matrix(
      data = 0,
      nrow = length(mods_name),
      ncol = length(mod_genes_name),
      byrow = FALSE,
      dimnames = NULL
    )
    row.names(mod_genes_table) = mods_name
    colnames(mod_genes_table) = mod_genes_name

    for (i in 1:length(rownames(mod_genes_table))) {
      tmp_col_mod = mods_name[i]
      tmp_mod_gene_json = module_gene_json[[i]]
      for (j in 1:length(colnames(mod_genes_table))) {
        tmp_row_gene = mod_genes_name[j]
        if (tmp_row_gene %in% tmp_mod_gene_json) {
          #mod_genes_table[i,j] = 1
          mod_genes_table[i, j] = genes_grad[tmp_row_gene, "grad"]
        }
      }

    }
    return(mod_genes_table)
  }

get_gene_df <-
  function(times,
           accurate,
           max_mgene_len,
           system_rightside,
           mod_genes_name,
           mod_genes_table) {
    #system_rightside = as.matrix(unlist(module_sim_value[1,]))
    gene_df = data.frame(rep(0, ncol(mod_genes_table)))
    gene_df = data.frame(colnames(mod_genes_table))

    #rank_sol = ncol(mod_genes_table)-nrow(mod_genes_table)
    for (i in 1:times) {
      if (i %% 1000 == 0) {
        print(i)
      }
      # add_random = runif(1,min=0.001,max=1.2)
      # bls = runif(length(mod_genes_name), min=0.00001, max=0.02)
      # bus = bls*1.2+add_random
      #add_random = runif(1,min=1.4999,max=1.5)
      add_random = runif(1,
                         min = 1.4999 * max_mgene_len,
                         max = 1.5 * max_mgene_len)
      bls = runif(length(mod_genes_name), min = 0.000001, max = 0.000002)
      bus = bls * 1.2 + add_random

      bvls_solution = bvls(mod_genes_table,
                           system_rightside,
                           bl = bls,
                           bu = bus)
      #sol = round(bvls_solution$x,6)
      sol = bvls_solution$x

      solution_bias = abs((mod_genes_table %*% sol) - system_rightside)
      if ((max(solution_bias) < accurate)) {
        gene_df = cbind(gene_df, sol)
      }
    }

    rownames(gene_df) = gene_df[, 1]
    gene_df = gene_df[, -1]
    print(length(gene_df))

    sol_correlation <- cor(gene_df)
    if (ncol(gene_df) > 0) {
      sample_nums = str_pad(1:ncol(gene_df), floor(log10(ncol(gene_df))) + 1, side =
                              "left", "0")
      colnames(gene_df) <- paste("sample_", sample_nums, sep = "")
    }
    return(gene_df)
  }


############################### gene_name #####################################
get_mod_genes_name <- function(module_gene_json) {
  # module_gene_json <- fromJSON(file = paste0("E:/scPlus_universal/simulate/flux_generater/org_input/",
  #                                            input_filename, "_modules_genes.json"))
  mod_genes_name = c()
  for (mod_gene in module_gene_json) {
    mod_genes_name = append(mod_genes_name, mod_gene)
  }

  mod_genes_name = unique(mod_genes_name)

  return(mod_genes_name)
}

############################### gene grad para #####################################
get_genes_grad <- function(mod_genes_name) {
  genes_grad = as.data.frame(mod_genes_name)
  rownames(genes_grad) = mod_genes_name
  colnames(genes_grad) = c("gene_name")
  grads = runif(length(mod_genes_name), min = 0.5, max = 1.5)
  #grads = runif(length(mod_genes_name),min=0.9,max=1.9)
  genes_grad[, "grad"] = grads
  return(genes_grad)
}


############################### generate #####################################
gene_expression_sim_main <-
  function(module_gene_json_path,
           module_sim_flux_path,
           times,
           accurate,
           csv_rowno_in_python = 0) {
    # setwd("E:/scPlus_universal/simulate/flux_generater/")

    # input_filename = "GSLsim1"
    # csv_rowno_in_python = 0
    # times = 80000

    # module_gene_json <- fromJSON(file = paste0("E:/scPlus_universal/simulate/flux_generater/org_input/",
    #                                            input_filename, "_modules_genes.json"))
    # module_sim_value <- read_csv(paste0("E:/scPlus_universal/simulate/flux_generater/simresult_flux/",
    #                                     input_filename, "_cmMat_sim.csv"), col_names = TRUE)

    module_gene_json <- fromJSON(file = module_gene_json_path)
    module_sim_value <-
      read_csv(module_sim_flux_path, col_names = TRUE)

    mod_genes_name = get_mod_genes_name(module_gene_json)
    genes_grad = get_genes_grad(mod_genes_name)
    if (ncol(genes_grad) > 0) {
      write.csv(
        genes_grad,
        paste0(
          "simresult_gene_grad/",
          input_filename,
          "_gene_grad_sim_row_",
          csv_rowno_in_python,
          ".csv"
        ),
        row.names = TRUE
      )
    }

    mgene_len = c()
    for (i in 1:length(module_gene_json)) {
      mgene_len = append(mgene_len, length(module_gene_json[[i]]))
    }
    max_mgene_len = max(mgene_len)

    ##
    mod_genes_table = get_mod_genes_table(module_gene_json, mod_genes_name, genes_grad)

    #m_gene_num = rowSums(as.matrix(mod_genes_table)==1)
    m_gene_num = rowSums(as.matrix(mod_genes_table) > 0)

    ##sim gene expression
    row_no = as.character(csv_rowno_in_python + 1)
    system_rightside = as.matrix(unlist(module_sim_value[row_no, ]))
    system_rightside_multi_mean = system_rightside

    for (i in 1:nrow(system_rightside_multi_mean)) {
      sys_module = names(system_rightside_multi_mean[i, 1])
      system_rightside_multi_mean[i, 1] = m_gene_num[sys_module] * system_rightside_multi_mean[i, 1]
    }

    gene_df = get_gene_df(
      times,
      accurate,
      max_mgene_len,
      system_rightside_multi_mean,
      mod_genes_name,
      mod_genes_table
    )
    # if (ncol(gene_df) > 0) {
    #   write.csv(gene_df, paste0("simresult_gene_expression/", input_filename, "_gene_expr_sim_row_", csv_rowno_in_python, ".csv"),row.names=TRUE)
    # }
    return(gene_df)
  }



############################### generate #####################################


#######################################################################################
#' sim_metabolism_flux function
#'
#' Generate Simulated Metabolic Flux Data.
#' @keywords sim_metabolism_flux
#' @param
#' sim_metabolism_flux:
#' Default:
#' precision=0.00000001
#' times=500
#' metabolism_network_path
#' @export
#' @examples
#' sim_metabolism_flux(precision=0.00000001, times=500, metabolism_network_path)
#'
#'
sim_metabolism_flux <-
  function(precision = 0.00000001,
           times = 500,
           metabolism_network_path) {
    # setwd("E:/scPlus_universal/simulate/flux_generater/")
    #csvpath = "E:/scPlus_universal/simulate/flux_generater/org_input/"
    #csvnames = c("GSL1_cmMat.csv", "GSL2_cmMat.csv", "GSLsim1_cmMat.csv")

    # set.seed(5)
    flux_df_sim = flux_sim_main(precision, times, metabolism_network_path)
    return(flux_df_sim)
  }



#######################################################################################
#' sim_metabolism_gene_expression function
#'
#' Generate Simulated Metabolic Flux Data.
#' @keywords sim_metabolism_gene_expression
#' @param
#' sim_metabolism_gene_expression:
#' Default:
#' module_gene_json_path
#' module_sim_flux_path
#' times=900
#' precision=0.00000000000001
#' metabolism_network_path
#' @export
#' @examples
#' sim_metabolism_gene_expression(module_gene_json_path, module_sim_flux_path, times=900, precision=0.00000000000001)
#'
#'
sim_metabolism_gene_expression <-
  function(module_gene_json_path,
           module_sim_flux_path,
           times = 900,
           precision = 0.00000000000001) {
    library(readr)
    library(bvls)
    library("rjson")
    library("stringr")


    # set.seed(5)
    # gene_df_sim = gene_expression_sim_main("GSL1", 0, 900, 0.00000000000001)
    gene_df_sim = gene_expression_sim_main(module_gene_json_path,
                                           module_sim_flux_path,
                                           times,
                                           precision,
                                           0)

    return(gene_df_sim)
  }
