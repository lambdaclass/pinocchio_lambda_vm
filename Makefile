test:
	cargo test

test_arkworks_adapter:
	cd arkworks_adapter && cargo test

docker-shell:
	docker build -t rust-curves .
	docker run -it rust-curves bash

nix-shell:
	nix-shell
