#!/bin/sh

head -6 new.txt | awk '{print $1}' > a
grep -w -f a *.name> b
