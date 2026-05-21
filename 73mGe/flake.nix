{
  description = "ROOT Analysis Development Environment";
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    agenix = {
      url = "github:ryantm/agenix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    utils = {
      url = "/home/e-work/Analysis-Utilities";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      agenix,
      utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
            cudaCapabilities = [ "12.0" ];
            cudaForwardCompat = false;
          };
        };
        analysis-utils = utils.packages.${system}.cuda;
        root = utils.packages.${system}.rootCuda;
        agenixPkg = agenix.packages.${system}.default;
      in
      {
        devShells.default = pkgs.mkShell {
          nativeBuildInputs = with pkgs; [
            pkg-config
            gnumake
            clang-tools
          ];
          buildInputs = [
            analysis-utils
            root
            agenixPkg
          ];
          shellHook = ''
            echo "Analysis-Utilities version: ${analysis-utils.version} (CUDA)"
            export NIX_CFLAGS_COMPILE="-DAU_ROOFIT_BACKEND_CUDA=1''${NIX_CFLAGS_COMPILE:+ $NIX_CFLAGS_COMPILE}"
            export CPLUS_INCLUDE_PATH="$PWD/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
            export ROOT_INCLUDE_PATH="$PWD/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
            export LD_LIBRARY_PATH="$PWD/lib:/run/opengl-driver/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
            alias clean-aclic='rm -f *_C.so *_C.d *_C_ACLiC_dict_rdict.pcm *_cpp.so *_cpp.d *_cpp_ACLiC_dict_rdict.pcm *_cxx.so *_cxx.d *_cxx_ACLiC_dict_rdict.pcm AutoDict_*'

            flake_root="$PWD"
            mkdir -p "$flake_root/include" "$flake_root/src" "$flake_root/macros"
            (
              cd "$flake_root/secrets"
              ${agenixPkg}/bin/agenix -d BEF.hpp.age -i "$HOME/.ssh/id_ed25519" > "$flake_root/include/BEF.hpp"
              ${agenixPkg}/bin/agenix -d BEF.cpp.age -i "$HOME/.ssh/id_ed25519" > "$flake_root/src/BEF.cpp"
              ${agenixPkg}/bin/agenix -d ConvertBEF.cpp.age -i "$HOME/.ssh/id_ed25519" > "$flake_root/macros/ConvertBEF.cpp"
            )
            chmod 644 "$flake_root/include/BEF.hpp" "$flake_root/src/BEF.cpp" "$flake_root/macros/ConvertBEF.cpp"

            cd macros
          '';
        };
      }
    );
}
