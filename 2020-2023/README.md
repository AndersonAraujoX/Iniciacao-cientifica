# Estados Quânticos de Dois Elétrons em Campo Magnético

Este repositório contém os códigos e materiais de apoio para o projeto de Iniciação Científica e Trabalho de Conclusão de Curso (TCC) intitulado "Cálculo e análise dos potenciais de interações dos elétrons em níveis de Landau num sistema bidimensional em campo magnético". O trabalho foi desenvolvido no Instituto de Física de São Carlos (IFSC) da Universidade de São Paulo (USP).

**Autor:** Anderson Araujo de Oliveira  
**Orientador:** Prof. Dr. Guo-Qiang Hai

## Resumo do Projeto

Neste projeto, estudamos os estados quânticos de dois elétrons em um sistema bidimensional, sujeitos a um campo magnético externo perpendicular e com interação Coulombiana entre eles. A equação de Schrödinger para este problema foi resolvida utilizando tanto métodos numéricos quanto soluções analíticas particulares (conhecidas como "quase exatas"). O objetivo é obter os espectros de energia, as funções de onda e comparar os resultados das diferentes abordagens para valores específicos do campo magnético.

## Estrutura do Repositório

  * **/relatorios**: Contém os relatórios do projeto.
      * `Relatório_PUB_2_corrigido.pdf`: Relatório final da Iniciação Científica (Agosto/2023).
      * `TCC_final.pdf`: Monografia do Trabalho de Conclusão de Curso (2022).
  * **/codigo\_fonte**: Contém os códigos-fonte utilizados nas simulações.
      * **Fortran** (`.f90`):
          * `nume.f90`: Implementação do método numérico para diagonalização da matriz Hamiltoniana.
          * `eigen.f90`: Sub-rotina para a diagonalização de matrizes (utiliza a biblioteca LAPACK).
          * `icq1.f90` e `icq2.f90`: Códigos para cálculo das soluções analíticas (quase exatas).
      * **Python/Jupyter** (`.ipynb`):
          * `ic.ipynb`: Notebook com cálculos e análises numéricas.
          * `plot_ic_2.ipynb` e `plot_ic_3.ipynb`: Notebooks para a geração dos gráficos e visualização dos resultados.

## Metodologia

A abordagem para resolver o problema foi dividida em duas frentes principais:

1.  **Solução Numérica:**

      * A Hamiltoniana do sistema é construída como uma matriz em uma base de estados de um elétron (níveis de Landau).
      * A matriz é diagonalizada numericamente para encontrar os autoestados (funções de onda) e autovalores (energias). O código `nume.f90` implementa essa lógica, utilizando a sub-rotina `eigen.f90` que por sua vez utiliza as rotinas `DSYEV` e `DSYEVD` da biblioteca LAPACK.

2.  **Solução Analítica (Quase Exata):**

      * Para valores particulares do campo magnético, existem soluções analíticas que podem ser expressas como polinômios finitos.
      * Os códigos `icq1.f90` e `icq2.f90` calculam essas soluções específicas, que servem como uma excelente forma de validar os resultados numéricos.

Os notebooks em Python foram utilizados para processar os dados gerados pelos códigos em Fortran, comparar os resultados dos dois métodos e gerar os gráficos presentes nos relatórios.

## Como Executar os Códigos

### Pré-requisitos

  * **Fortran:** Um compilador Fortran (como o `gfortran`) e a biblioteca `LAPACK`.
      * Para compilar: `gfortran -o executavel nome_do_arquivo.f90 -llapack`
  * **Python:** Python 3 com as bibliotecas `numpy`, `scipy`, `matplotlib` e `jupyter`.
      * Para executar os notebooks: `jupyter notebook`

### Passos

1.  Clone o repositório:
    ```bash
    git clone [URL-DO-SEU-REPOSITORIO]
    ```
2.  Compile e execute os códigos em Fortran para gerar os dados numéricos e analíticos.
3.  Abra os notebooks Jupyter (`.ipynb`) para visualizar e analisar os resultados.

## Resultados

Os resultados, incluindo os espectros de energia, as funções de onda e a comparação detalhada entre os métodos numérico e analítico, podem ser encontrados nos documentos da pasta `/relatorios`. Os gráficos demonstram uma excelente concordância entre as duas abordagens nos casos onde a solução analítica é aplicável.

## Referências Principais

  * Taut, M. Two electrons in a homogeneous magnetic field: particular analytical solutions. *Journal of Physics A: Mathematical and General*, 27(3):1045, 1994.
  * Turbiner, A. V. Quasi-exactly-solvable problems and sl(2) algebra. *Communications in Mathematical Physics*, 118(3):467-474, 1988.

*(Para uma lista completa, consulte a seção de referências nos relatórios.)*
