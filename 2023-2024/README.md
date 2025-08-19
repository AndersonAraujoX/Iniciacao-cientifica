# Tunelamento Dissipativo Não-Local

**Relatório de Iniciação Científica (IC) apresentado ao Programa Unificado de Bolsas de Estudo para Estudantes de Graduação da Universidade de São Paulo.**

- **Bolsista:** Anderson Araujo de Oliveira
- **Orientador:** Prof. Dr. Miled Hassan Youssef Moussa
- **Instituição:** Instituto de Física de São Carlos (IFSC), Universidade de São Paulo

## Resumo

Este projeto de pesquisa foca no estudo do processo de tunelamento dissipativo não local, que é um mecanismo para a transferência de estados em redes lineares de osciladores dissipativos. Durante o período de seis meses da bolsa, o aluno se dedicou ao estudo da integral de trajetória e do modelo de Caldeira-Leggett. O trabalho culminou na dedução da equação de movimento para o referido modelo.

## Objetivos do Projeto

Os principais objetivos desta pesquisa foram:

* Analisar o tunelamento não-local de uma partícula entre dois poços de potencial distantes, conectados por um terceiro poço.
* Desenvolver um tratamento mais realista para o tunelamento de uma partícula entre poços distantes, considerando as posições da partícula.
* Adaptar o método de redes de osciladores para os modelos de Caldeira-Leggett e Spin-Boson.

## Metodologia

A metodologia empregada neste estudo baseou-se na formulação de integrais de trajetória de Feynman para descrever a evolução de sistemas quânticos. Foi utilizado o conceito de propagador (função de Green dependente do tempo) para calcular a evolução do sistema entre dois instantes. A abordagem foi desenvolvida tanto no espaço de fase quanto no espaço de configurações.

## Resultados

O estudo se concentrou no modelo de Caldeira-Leggett, que descreve um sistema composto por um sistema principal, um reservatório (formado por um conjunto de osciladores) e a interação entre eles. A partir da Lagrangiana do sistema, foi derivada a equação de Langevin.

Utilizando a abordagem do superpropagador, que descreve a evolução temporal do operador de densidade reduzido do sistema, e o funcional de influência de Feynman-Vernon, foi possível obter a equação de movimento para o operador de densidade do sistema no limite de altas temperaturas (região semiclássica).

A equação de movimento final para a matriz de densidade reduzida ($\rho_s$) na região semiclássica é dada por:

$$\frac{\partial\rho_{s}}{\partial t}=-\frac{\hbar}{2Mi}\frac{\partial^{2}\rho_{s}}{\partial x^{2}}+\frac{\hbar}{2Mi}\frac{\partial^{2}\rho_{s}}{\partial y^{2}}-\gamma(x-y)[\frac{\partial\rho_{s}}{\partial x}-\frac{\partial\rho_{s}}{\partial y}]+\frac{1}{i\hbar}[\tilde{V}(x)-\tilde{V}(y)]\rho_{s}-\frac{2MVk_{B}T}{\hbar^{2}}(x-y)^{2}p_{s}$$

## Conclusão

O bolsista logrou êxito em seus estudos sobre a integral de trajetória e o modelo de Caldeira-Leggett, conseguindo reproduzir os cálculos e deduzir a equação de movimento para o modelo em questão. O trabalho desenvolvido servirá de base para a continuação da pesquisa sobre tunelamento não-local em nível de mestrado.

## Referências

O relatório cita trabalhos de referência nas áreas de transferência de estados quânticos, informação quântica e sistemas quânticos dissipativos, incluindo:

* **Caldeira, A. O., & Leggett, A. J. (1983).** Quantum tunnelling in a dissipative system. *Annals of physics*, 149(2), 374-456.
* **Neto, G. D. M., de Ponte, M. A., & Moussa, M. H. Y. (2012).** Nonlocal dissipative tunneling for high-fidelity quantum-state transfer between distant parties. *Physical Review A*, 85(5), 052303.
