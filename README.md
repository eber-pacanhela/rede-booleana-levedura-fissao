# Implementação em R de uma rede booleana para a predição da sequência do ciclo celular da levedura de fissão

Este repositório contém uma implementação, em linguagem R, da rede booleana proposta por Davidich e Bornholdt (2008) para representar a dinâmica regulatória do ciclo celular da levedura de fissão (*Schizosaccharomyces pombe*).

O projeto foi desenvolvido como parte do relatório técnico da disciplina de **Biologia de Sistemas**, vinculada ao **Programa de Pós-Graduação Associado em Bioinformática da UFPR/UTFPR**.

## Objetivo

O objetivo do projeto é reconstruir computacionalmente o modelo booleano descrito no artigo de referência, permitindo:

- representar os principais componentes e as interações regulatórias do ciclo celular da levedura de fissão;
- implementar a regra de atualização dos estados dos nós;
- converter automaticamente a dinâmica do modelo em funções compatíveis com o pacote `BoolNet`;
- reproduzir a sequência temporal do ciclo celular;
- identificar exaustivamente os atratores e suas respectivas bacias de atração;
- verificar a correspondência entre a implementação direta da regra de atualização e as transições calculadas pelo `BoolNet`.

## Artigo de referência

> Davidich, M. I.; Bornholdt, S. Boolean Network Model Predicts Cell Cycle Sequence of Fission Yeast. *PLOS ONE*, 3, 2008.
> DOI: [10.1371/journal.pone.0001672](https://doi.org/10.1371/journal.pone.0001672)

## Componentes da rede

A rede é formada por dez nós:

- `Start`;
- `SK`;
- `Cdc2_Cdc13`;
- `Ste9`;
- `Rum1`;
- `Slp1`;
- `Cdc2_Cdc13_estrela`;
- `Wee1_Mik1`;
- `Cdc25`;
- `PP`.

Cada nó assume um estado binário:

- `0`: componente ausente ou inativo;
- `1`: componente presente ou ativo.

As interações são armazenadas em uma matriz quadrada de dimensão 10 × 10, na qual:

- `+1` representa uma interação de ativação;
- `-1` representa uma interação de inibição;
- `0` representa ausência de interação.

Os nós `Start`, `SK`, `Slp1` e `PP` possuem laços negativos de autodegradação.

## Regra de atualização

Para cada nó, calcula-se a soma ponderada das entradas regulatórias. O estado seguinte é determinado pela comparação entre essa soma e o limiar do nó:

- o nó é ativado quando a soma é superior ao limiar;
- o nó é inativado quando a soma é inferior ao limiar;
- o estado anterior é mantido quando a soma é igual ao limiar.

O limiar é igual a zero para a maioria dos componentes. Foram utilizados os seguintes valores especiais:

- `Cdc2_Cdc13`: `-0,5`;
- `Cdc2_Cdc13_estrela`: `0,5`.

A atualização é síncrona, isto é, todos os nós são atualizados simultaneamente em cada passo de tempo.

## Funcionalidades do programa

O script realiza as seguintes etapas:

1. instala e carrega o pacote `BoolNet`, caso necessário;
2. define os dez nós da rede;
3. constrói a matriz de interações regulatórias;
4. aplica os limiares específicos de ativação;
5. implementa a regra de atualização dos estados;
6. gera automaticamente as funções booleanas de cada nó;
7. grava as regras no arquivo `regras_levedura_fissao.txt`;
8. carrega a rede no `BoolNet`;
9. valida todas as 1.024 transições possíveis do espaço de estados;
10. simula a trajetória biológica a partir da condição inicial definida;
11. executa uma busca exaustiva pelos atratores;
12. apresenta os atratores e os tamanhos de suas bacias no console.

## Requisitos

A implementação foi desenvolvida com:

- R 4.6.0;
- pacote `BoolNet` 2.1.9.

O próprio script instala o pacote `BoolNet` automaticamente quando ele não está disponível no ambiente.

## Execução

Clone o repositório:

```bash
git clone https://github.com/eber-pacanhela/rede-booleana-levedura-fissao.git
```

Abra o arquivo `codigo.r` no RStudio e execute-o integralmente.

## Condição inicial

A simulação parte do seguinte estado biológico:

```text
Start = 1
SK = 0
Cdc2_Cdc13 = 0
Ste9 = 1
Rum1 = 1
Slp1 = 0
Cdc2_Cdc13_estrela = 0
Wee1_Mik1 = 1
Cdc25 = 0
PP = 0
```

Portanto, inicialmente estão ativos `Start`, `Ste9`, `Rum1` e `Wee1_Mik1`, enquanto os demais componentes estão inativos.

## Resultados esperados

A execução do programa deve:

- confirmar a correspondência das 1.024 transições calculadas diretamente e pelo `BoolNet`;
- reproduzir uma trajetória temporal de dez estados, com retorno ao estado estacionário G1;
- identificar 13 atratores;
- obter bacias de atração com tamanhos 762, 208, 18, 18 e nove bacias de tamanho 2;
- confirmar que a soma dos tamanhos das bacias corresponde aos 1.024 estados possíveis da rede.

O ponto fixo dominante corresponde ao estado estacionário G1, no qual `Ste9`, `Rum1` e `Wee1_Mik1` permanecem ativos.

## Arquivos

```text
.
├── codigo.r
└── README.md
```

## Reprodutibilidade

A dinâmica implementada é determinística e utiliza atualização síncrona. Assim, para a mesma versão do código e das dependências, a execução deve produzir a mesma trajetória, os mesmos atratores e os mesmos tamanhos de bacia.

## Contexto acadêmico

**Instituições**

- Universidade Federal do Paraná - UFPR;
- Universidade Tecnológica Federal do Paraná - UTFPR;
- Programa de Pós-Graduação Associado em Bioinformática.

**Disciplina**

- Biologia de Sistemas.

**Docentes responsáveis**

- Profa. Dra. Deisy Gysi;
- Prof. Dr. Fabrício Lopes;
- Prof. Dr. Mauro Castro.

## Autor

**Eber Fabiano Pacanhela**

Cornélio Procópio, 2026.
