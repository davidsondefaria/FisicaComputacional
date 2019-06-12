#!/bin/bash
./proj >> dados.dat
cat dados.dat|awk '{print $2 " "$3 " "$4}' >> energia.dat
cat dados.dat|awk '{print $5 " "$6 " "$7 " " $8}' >> posicoes.dat
