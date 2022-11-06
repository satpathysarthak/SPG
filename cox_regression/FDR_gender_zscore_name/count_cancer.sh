#!/bin/sh

head -8 new.txt | awk '{print $1}' > a
grep -w -f a *.name> b
