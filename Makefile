test:
	cargo test

docker-shell:
	docker build -t cuda118 .
	docker run --rm -it --mount src=$(PWD),dst=/cuda-code,type=bind cuda118

nix-shell:
	nix-shell
