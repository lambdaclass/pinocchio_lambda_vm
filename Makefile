run:
	cargo run --release

docker-shell:
	docker build -t cuda118 .
	docker run --rm -it --mount src=$(PWD),dst=/cuda-code,type=bind cuda118

docker-cuda-shell:
	docker build -t cuda118 .
	docker run --rm --gpus all -it --mount src=$(PWD),dst=/cuda-code,type=bind cuda118

nix-shell:
	nix-shell

nix-cuda-shell:
	NIXPKGS_ALLOW_UNFREE=1 nix-shell --impure cuda-shell.nix
