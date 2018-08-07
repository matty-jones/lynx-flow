#!/bin/bash
rsync -r --exclude="*.py" --exclude="*.o*" --exclude="*.e*" --exclude="*.OU" --exclude ="*.ER" --exclude="*.sh" --exclude="*.md" --exclude="*.txt" --exclude="*pycache*" --update --stats --progress anderso4@bw.ncsa.illinois.edu://u/eot/anderso4/projects/rhaco-flow/* ./
