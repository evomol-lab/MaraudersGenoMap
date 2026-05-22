PYTHON = python3
GUI = marauders.py
BIN = /usr/local/bin/marauders
IMAGE = evomol/marauders-genomap:latest
SIF = marauders-genomap.sif
APP_DIR = /app

run:
	@$(PYTHON) $(GUI)

docker-build:
	docker build --no-cache -t $(IMAGE) .

docker-build-cache:
	docker build -t $(IMAGE) .

docker-run:
	docker run --rm -it \
		--network host \
		--user $$(id -u):$$(id -g) \
		-e DISPLAY=$$DISPLAY \
		-e XAUTHORITY=/tmp/.Xauthority \
		-e HOME=/tmp \
		-e QT_X11_NO_MITSHM=1 \
		-v /tmp/.X11-unix:/tmp/.X11-unix:rw \
		-v "$${XAUTHORITY:-$$HOME/.Xauthority}":/tmp/.Xauthority:ro \
		-v "$(PWD)":$(APP_DIR) \
		-w $(APP_DIR) \
		$(IMAGE)

docker-run-ssh:
	$(MAKE) docker-run

install:
	@echo "#!/bin/bash" > marauders_cmd
	@echo "cd $(PWD) && $(PYTHON) $(GUI)" >> marauders_cmd
	@chmod +x marauders_cmd
	@sudo mv marauders_cmd $(BIN)
	@echo "✅ Agora você pode abrir o programa digitando apenas: marauders"

uninstall:
	@sudo rm -f $(BIN)
