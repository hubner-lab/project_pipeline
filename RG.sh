#!/bin/bash
sm=$(echo $1 | cut -d'_' -f1-2)
id="ID"
echo "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA"
