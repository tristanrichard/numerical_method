# Nom de l'exécutable
EXEC = bin\run.exe

# Dossiers
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Fichiers source
SRC = $(wildcard $(SRC_DIR)/*.c)

# Fichiers objets
OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# Compilateur et options
CC = gcc
CFLAGS = -Wall -Wextra -Iheader

# Règle pour la compilation de l'exécutable
$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $(EXEC)

# Règle pour la compilation des fichiers .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# Règle pour exécuter le programme avec des arguments (recompilation à chaque fois)
.PHONY: run
run: clean $(EXEC)
	./$(EXEC) $(ARGS)

# Règle pour nettoyer les fichiers objets et l'exécutable
clean:
	del /f $(OBJ_DIR)\*.o $(EXEC)

.PHONY: clean