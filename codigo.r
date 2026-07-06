# ==============================================================
#           Implementação em R de uma rede booleana         
#        para a predição da sequência do ciclo celular 
#                   da levedura de fissão
#
# Implementação desenvolvida com base na descrição das 
# interações regulatórias, na dinâmica de atualização, 
# na condição inicial biológica e nos resultados apresentados 
# por Davidich e Bornholdt (2008).
#
# Artigo de referência:
# Boolean Network Model Predicts Cell Cycle Sequence of Fission Yeast
# DOI: 10.1371/journal.pone.0001672 
# ==============================================================


# --------------------------------------------------------------
# 1. Instalação e carregamento do pacote
# --------------------------------------------------------------

if (!requireNamespace("BoolNet", quietly = TRUE)) {
  install.packages("BoolNet")
}

library(BoolNet)


# --------------------------------------------------------------
# 2. Nós da rede
# --------------------------------------------------------------

nos <- c(
  "Start",
  "SK",
  "Cdc2_Cdc13",
  "Ste9",
  "Rum1",
  "Slp1",
  "Cdc2_Cdc13_estrela",
  "Wee1_Mik1",
  "Cdc25",
  "PP"
)

numero_nos <- length(nos)


# --------------------------------------------------------------
# 3. Matriz de interações regulatórias
#
# Linhas: nós que recebem a interação
# Colunas: nós que enviam a interação
#
#  1 = ativação
# -1 = inibição
#  0 = ausência de interação
# --------------------------------------------------------------

matriz_interacoes <- matrix(
  0L,
  nrow = numero_nos,
  ncol = numero_nos,
  dimnames = list(
    alvo = nos,
    origem = nos
  )
)


adicionar_interacao <- function(origem, alvo, peso) {

  if (!(origem %in% nos)) {
    stop("Nó de origem inexistente: ", origem)
  }

  if (!(alvo %in% nos)) {
    stop("Nó de destino inexistente: ", alvo)
  }

  if (!(peso %in% c(-1L, 1L))) {
    stop("O peso deve ser -1 ou 1.")
  }

  matriz_interacoes[alvo, origem] <<- as.integer(peso)
}


# --------------------------------------------------------------
# 3.1 Interações de ativação
# --------------------------------------------------------------

# Start ativa as quinases iniciais.
adicionar_interacao("Start", "SK", 1L)

# Cdc2/Cdc13 ativa Cdc25.
adicionar_interacao("Cdc2_Cdc13", "Cdc25", 1L)

# Cdc2/Cdc13* ativa Slp1.
adicionar_interacao("Cdc2_Cdc13_estrela", "Slp1", 1L)

# Slp1 ativa PP.
adicionar_interacao("Slp1", "PP", 1L)

# PP ativa Ste9, Rum1 e Wee1/Mik1.
adicionar_interacao("PP", "Ste9", 1L)
adicionar_interacao("PP", "Rum1", 1L)
adicionar_interacao("PP", "Wee1_Mik1", 1L)

# Cdc25 ativa a forma altamente ativa de Cdc2/Cdc13.
adicionar_interacao("Cdc25", "Cdc2_Cdc13_estrela", 1L)


# --------------------------------------------------------------
# 3.2 Interações de inibição
# --------------------------------------------------------------

# SK inibe Ste9 e Rum1.
adicionar_interacao("SK", "Ste9", -1L)
adicionar_interacao("SK", "Rum1", -1L)

# Wee1/Mik1 inibe Cdc2/Cdc13*.
adicionar_interacao("Wee1_Mik1", "Cdc2_Cdc13_estrela", -1L)

# Rum1 inibe Cdc2/Cdc13 e Cdc2/Cdc13*.
adicionar_interacao("Rum1", "Cdc2_Cdc13", -1L)
adicionar_interacao("Rum1", "Cdc2_Cdc13_estrela", -1L)

# Cdc2/Cdc13 inibe Rum1.
adicionar_interacao("Cdc2_Cdc13", "Rum1", -1L)

# Ste9 inibe Cdc2/Cdc13 e Cdc2/Cdc13*.
adicionar_interacao("Ste9", "Cdc2_Cdc13", -1L)
adicionar_interacao("Ste9", "Cdc2_Cdc13_estrela", -1L)

# Slp1 inibe Cdc2/Cdc13 e Cdc2/Cdc13*.
adicionar_interacao("Slp1", "Cdc2_Cdc13", -1L)
adicionar_interacao("Slp1", "Cdc2_Cdc13_estrela", -1L)

# Cdc2/Cdc13 inibe Ste9 e Wee1/Mik1.
adicionar_interacao("Cdc2_Cdc13", "Ste9", -1L)
adicionar_interacao("Cdc2_Cdc13", "Wee1_Mik1", -1L)

# PP inibe Cdc25.
adicionar_interacao("PP", "Cdc25", -1L)

# Cdc2/Cdc13* inibe Rum1 e Ste9.
adicionar_interacao("Cdc2_Cdc13_estrela", "Rum1", -1L)
adicionar_interacao("Cdc2_Cdc13_estrela", "Ste9", -1L)


# --------------------------------------------------------------
# 3.3 Autodegradação
#
# Os nós Start, SK, Slp1 e PP possuem laços negativos.
# --------------------------------------------------------------

adicionar_interacao("Start", "Start", -1L)
adicionar_interacao("SK", "SK", -1L)
adicionar_interacao("Slp1", "Slp1", -1L)
adicionar_interacao("PP", "PP", -1L)


# --------------------------------------------------------------
# 4. Limiares de ativação
#
# Como as somas das entradas são inteiras, os valores -0,5 e
# 0,5 representam os limiares especiais necessários para
# reproduzir a dinâmica apresentada no artigo.
# --------------------------------------------------------------

limiares <- setNames(
  rep(0, numero_nos),
  nos
)

limiares["Cdc2_Cdc13"] <- -0.5
limiares["Cdc2_Cdc13_estrela"] <- 0.5


# --------------------------------------------------------------
# 5. Implementação da regra de atualização dos estados
#
# soma_i(t) = soma_j a_ij * S_j(t)
#
# S_i(t+1) =
#   1,      se soma_i(t) > limiar_i
#   0,      se soma_i(t) < limiar_i
#   S_i(t), se soma_i(t) = limiar_i
# --------------------------------------------------------------

aplicar_regra_atualizacao <- function(estado) {

  if (length(estado) != numero_nos) {
    stop(
      "O estado deve possuir ",
      numero_nos,
      " valores."
    )
  }

  if (is.null(names(estado))) {
    names(estado) <- nos
  }

  estado <- as.integer(estado[nos])
  names(estado) <- nos

  if (any(!estado %in% c(0L, 1L))) {
    stop("Todos os estados devem ser 0 ou 1.")
  }

  somas <- drop(matriz_interacoes %*% estado)

  proximo_estado <- estado

  proximo_estado[somas > limiares] <- 1L
  proximo_estado[somas < limiares] <- 0L

  # Quando soma == limiar, o valor anterior é mantido.

  names(proximo_estado) <- nos

  return(proximo_estado)
}


# --------------------------------------------------------------
# 6. Conversão automática da regra de atualização dos estados 
# para regras do BoolNet
#
# Para cada nó são enumeradas todas as combinações possíveis
# de seus reguladores. As combinações que produzem estado 1
# são convertidas para uma expressão na forma normal disjuntiva.
# --------------------------------------------------------------

gerar_regra_boolnet <- function(alvo) {

  reguladores <- colnames(matriz_interacoes)[
    matriz_interacoes[alvo, ] != 0L
  ]

  # O próprio nó precisa ser considerado porque a regra de 
  # atualização dos estados mantém seu valor anterior quando a 
  # soma é igual ao limiar.
  entradas <- unique(c(reguladores, alvo))

  combinacoes <- expand.grid(
    rep(list(c(0L, 1L)), length(entradas)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  names(combinacoes) <- entradas

  saidas <- apply(
    combinacoes,
    1,
    function(linha) {

      estado <- setNames(
        integer(numero_nos),
        nos
      )

      estado[entradas] <- as.integer(linha)

      aplicar_regra_atualizacao(estado)[alvo]
    }
  )

  # Função constantemente falsa.
  if (all(saidas == 0L)) {
    return(
      paste0(
        "(",
        alvo,
        " & !",
        alvo,
        ")"
      )
    )
  }

  # Função constantemente verdadeira.
  if (all(saidas == 1L)) {
    return(
      paste0(
        "(",
        alvo,
        " | !",
        alvo,
        ")"
      )
    )
  }

  combinacoes_ativas <- combinacoes[
    saidas == 1L,
    ,
    drop = FALSE
  ]

  termos <- apply(
    combinacoes_ativas,
    1,
    function(linha) {

      literais <- ifelse(
        as.integer(linha) == 1L,
        names(linha),
        paste0("!", names(linha))
      )

      paste0(
        "(",
        paste(literais, collapse = " & "),
        ")"
      )
    }
  )

  paste(termos, collapse = " | ")
}


regras_booleanas <- vapply(
  nos,
  gerar_regra_boolnet,
  character(1)
)


# --------------------------------------------------------------
# 7. Geração do arquivo de regras
# --------------------------------------------------------------

arquivo_regras <- "regras_levedura_fissao.txt"

conteudo_arquivo <- c(
  "targets, functions",
  paste(
    nos,
    regras_booleanas,
    sep = ", "
  )
)

writeLines(
  conteudo_arquivo,
  arquivo_regras
)

cat(
  "\nArquivo de regras gerado:",
  arquivo_regras,
  "\n"
)


# --------------------------------------------------------------
# 8. Carregamento da rede pelo BoolNet
# --------------------------------------------------------------

rede <- BoolNet::loadNetwork(
  arquivo_regras
)

# O BoolNet identifica funções constantes como genes fixos.
# Entretanto, Start precisa participar dos 2^10 estados iniciais
# para que os tamanhos das bacias sejam os mesmos do artigo.
rede$fixed[] <- -1L

cat("\nRede carregada pelo BoolNet:\n")
print(rede)


# --------------------------------------------------------------
# 9. Validação das regras
#
# Compara as transições do BoolNet com a regra de atualização 
# dos estados para todos os 2^10 = 1.024 estados possíveis.
# --------------------------------------------------------------

todos_estados <- expand.grid(
  rep(list(c(0L, 1L)), numero_nos),
  KEEP.OUT.ATTRS = FALSE
)

todos_estados <- as.matrix(todos_estados)
colnames(todos_estados) <- nos

regras_validas <- all(
  apply(
    todos_estados,
    1,
    function(estado) {

      estado <- setNames(
        as.integer(estado),
        nos
      )

      resultado_equacao <- aplicar_regra_atualizacao(estado)

      resultado_boolnet <- BoolNet::stateTransition(
        rede,
        estado[rede$genes],
        type = "synchronous"
      )

      all(
        resultado_equacao[rede$genes] ==
          resultado_boolnet
      )
    }
  )
)

if (!regras_validas) {
  stop(
    "As regras do BoolNet não correspondem a regra de atualização."
  )
}

cat(
  "\nValidação concluída: as 1.024 transições",
  "correspondem a regra de atualização.\n"
)


# --------------------------------------------------------------
# 10. Condição inicial biológica
#
# Ativos:
#   Start, Ste9, Rum1 e Wee1/Mik1
#
# Inativos:
#   demais elementos
# --------------------------------------------------------------

estado_inicial <- c(
  Start = 1L,
  SK = 0L,
  Cdc2_Cdc13 = 0L,
  Ste9 = 1L,
  Rum1 = 1L,
  Slp1 = 0L,
  Cdc2_Cdc13_estrela = 0L,
  Wee1_Mik1 = 1L,
  Cdc25 = 0L,
  PP = 0L
)


# --------------------------------------------------------------
# 11. Simulação da sequência temporal
# --------------------------------------------------------------

numero_passos <- 10L

trajetoria <- matrix(
  0L,
  nrow = numero_passos,
  ncol = numero_nos,
  dimnames = list(
    NULL,
    nos
  )
)

trajetoria[1, ] <- estado_inicial[nos]

for (passo in 2:numero_passos) {

  trajetoria[passo, ] <- BoolNet::stateTransition(
    rede,
    trajetoria[passo - 1, ],
    type = "synchronous"
  )
}

fases <- c(
  "START",
  "G1",
  "G1/S",
  "G2",
  "G2",
  "G2/M",
  "G2/M",
  "M",
  "M",
  "G1"
)

tabela_temporal <- data.frame(
  Passo = seq_len(numero_passos),
  trajetoria,
  Fase = fases,
  check.names = FALSE
)

cat(
  "\n====================================================\n",
  "SEQUÊNCIA TEMPORAL DO CICLO CELULAR\n",
  "====================================================\n",
  sep = ""
)

print(
  tabela_temporal,
  row.names = FALSE
)


# --------------------------------------------------------------
# 12. Busca exaustiva dos atratores
#
# Todos os 1.024 estados são processados pelo BoolNet.
# --------------------------------------------------------------

atratores <- BoolNet::getAttractors(
  rede,
  type = "synchronous",
  method = "exhaustive",
  returnTable = TRUE
)

tamanhos_bacias <- vapply(
  atratores$attractors,
  function(atrator) {
    atrator$basinSize
  },
  numeric(1)
)

ordem_atratores <- order(
  tamanhos_bacias,
  decreasing = TRUE
)


# --------------------------------------------------------------
# 13. Construção da tabela de atratores e tamanhos de bacia
# --------------------------------------------------------------

linhas_tabela_atratores <- lapply(
  seq_along(ordem_atratores),
  function(numero_atrator) {

    indice_boolnet <- ordem_atratores[numero_atrator]

    estados_atrator <- BoolNet::getAttractorSequence(
      atratores,
      indice_boolnet
    )

    quantidade_estados <- nrow(estados_atrator)

    tipo_atrator <- if (
      quantidade_estados == 1L
    ) {
      "FP"
    } else {
      "LC"
    }

    data.frame(
      Atrator = rep(
        numero_atrator,
        quantidade_estados
      ),
      Tipo = rep(
        tipo_atrator,
        quantidade_estados
      ),
      Tamanho_bacia = rep(
        tamanhos_bacias[indice_boolnet],
        quantidade_estados
      ),
      estados_atrator,
      check.names = FALSE
    )
  }
)

tabela_atratores <- do.call(
  rbind,
  linhas_tabela_atratores
)

rownames(tabela_atratores) <- NULL

cat(
  "\n====================================================\n",
  "ATRATORES E TAMANHOS DAS BACIAS\n",
  "====================================================\n",
  sep = ""
)

print(
  tabela_atratores,
  row.names = FALSE
)