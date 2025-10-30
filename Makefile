
SHELL=/bin/bash

LUAVER=5.4
NAME=lambert93

CC=gcc

DBGFLAG ?= -Ofast -march=native -mtune=native -funroll-loops

LUAINCL=$(shell pkg-config --cflags lua$(LUAVER))
LUALIB=$(shell pkg-config --libs lua$(LUAVER))

CFLAGS= -std=gnu17 -D_GNU_SOURCE -D_REENTRANT $(DBGFLAG) -shared -fPIC -Wall -Wno-parentheses -Wno-switch -Wno-pointer-sign -Wno-trampolines -Wno-unused-result $(LUAINCL)
LIBS=-lm

SRC= $(wildcard *.c)

TARGET=$(NAME).so

.PHONY: clean dist install all

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $^ $(CFLAGS) $(LIBS) $(LUALIB) -o $@

rock:
	luarocks make $(NAME)-scm-1.rockspec

dist clean:
	rm -fR $(TARGET) semantic.cache* *.tmp *.tmp~ docs

