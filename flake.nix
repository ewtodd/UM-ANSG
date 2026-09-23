{
  description = "ANSG analysis projects: one dev shell per project";
  inputs = {
    utils.url = "github:ewtodd/Analysis-Utilities";
    nixpkgs.follows = "utils/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
    agenix = {
      url = "github:ryantm/agenix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      utils,
      agenix,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = {
            allowUnfree = true;
            cudaCapabilities = [ "8.9" ];
            cudaForwardCompat = false;
            problems.handlers = {
              torch.unsupported-cuda-version = "ignore";
            };
          };
          overlays = [
            (final: prev: {
              pythonPackagesExtensions = prev.pythonPackagesExtensions ++ [
                (python-final: python-prev: {
                  slicer = python-prev.slicer.overridePythonAttrs (old: {
                    doCheck = false;
                    nativeBuildInputs = (old.nativeBuildInputs or [ ]) ++ [
                      python-final.setuptools
                    ];
                  });
                  shap = python-prev.shap.overridePythonAttrs (old: {
                    dependencies = (old.dependencies or [ ]) ++ [
                      python-final.typing-extensions
                    ];
                  });
                  torch-bin =
                    (python-prev.torch-bin.override { cudaPackages = final.cudaPackages_13; })
                    .overridePythonAttrs
                      (old: {
                        pythonRelaxDeps = (old.pythonRelaxDeps or [ ]) ++ [
                          "setuptools"
                          "cuda-bindings"
                        ];
                      });
                })
              ];
            })
          ];
        };
        au = utils.packages.${system};
        agenixPkg = agenix.packages.${system}.default;

        tools = with pkgs; [
          pkg-config
          gnumake
          clang-tools
        ];

        # Every shell: locate the checkout, point the compilers at the
        # project's include/ and lib/, keep the ACLiC cleanup alias.
        projectHook = project: ''
          repo_root="$(git rev-parse --show-toplevel 2>/dev/null || echo "$PWD")"
          project_root="$repo_root/${project}"
          export CPLUS_INCLUDE_PATH="$project_root/include''${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"
          export ROOT_INCLUDE_PATH="$project_root/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
          export LD_LIBRARY_PATH="$project_root/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
          alias clean-aclic='rm -f *_C.so *_C.d *_C_ACLiC_dict_rdict.pcm *_cpp.so *_cpp.d *_cpp_ACLiC_dict_rdict.pcm *_cxx.so *_cxx.d *_cxx_ACLiC_dict_rdict.pcm AutoDict_*'
        '';

        analysisShell =
          {
            project,
            analysis-utils,
            root,
            extraInputs ? [ ],
            extraHook ? "",
            enterMacros ? true,
          }:
          pkgs.mkShell {
            nativeBuildInputs = tools;
            buildInputs = [
              analysis-utils
              root
              pkgs.bash
            ]
            ++ extraInputs;
            shellHook = ''
              echo "${project}: Analysis-Utilities ${analysis-utils.version}"
              ${projectHook project}
              ${extraHook}
              ${if enterMacros then ''cd "$project_root/macros"'' else ''cd "$project_root"''}
            '';
          };
      in
      {
        devShells = {
          # CUDA RooFit backend, ROOT built with CUDA, agenix-decrypted BEF
          # reader sources. Decrypts with the user's ssh key on entry.
          "73mGe" = analysisShell {
            project = "73mGe";
            analysis-utils = au.cuda;
            root = au.rootCuda;
            extraInputs = [ agenixPkg ];
            extraHook = ''
              export NIX_CFLAGS_COMPILE="-DAU_ROOFIT_BACKEND_CUDA=1''${NIX_CFLAGS_COMPILE:+ $NIX_CFLAGS_COMPILE}"
              export LD_LIBRARY_PATH="/run/opengl-driver/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
              mkdir -p "$project_root/include" "$project_root/src" "$project_root/macros"
              (
                cd "$project_root/secrets"
                ${agenixPkg}/bin/agenix -d BEF.hpp.age -i "$HOME/.ssh/id_ed25519" > "$project_root/include/BEF.hpp"
                ${agenixPkg}/bin/agenix -d BEF.cpp.age -i "$HOME/.ssh/id_ed25519" > "$project_root/src/BEF.cpp"
                ${agenixPkg}/bin/agenix -d ConvertBEF.cpp.age -i "$HOME/.ssh/id_ed25519" > "$project_root/macros/ConvertBEF.cpp"
              )
              chmod 644 "$project_root/include/BEF.hpp" "$project_root/src/BEF.cpp" "$project_root/macros/ConvertBEF.cpp"
            '';
          };

          "78mBr" = analysisShell {
            project = "78mBr";
            analysis-utils = au.default;
            root = pkgs.root;
          };

          "YAG-PreProcessing" = analysisShell {
            project = "YAG-PreProcessing";
            analysis-utils = au.default;
            root = pkgs.root;
            enterMacros = false;
            extraInputs = [
              (pkgs.python3.withPackages (ps: [
                ps.numpy
                ps.pandas
                ps.h5py
                au.pythonPackage
              ]))
            ];
          };

          "YAP-PSD" = analysisShell {
            project = "YAP-PSD";
            analysis-utils = au.default;
            root = pkgs.root;
            enterMacros = false;
            extraInputs = [
              (pkgs.python3.withPackages (
                ps: with ps; [
                  numpy
                  pandas
                  scikit-learn
                  xgboost
                  shap
                  packaging
                  torch-bin
                  au.pythonPackage
                ]
              ))
            ];
          };

          default = pkgs.mkShell {
            shellHook = ''
              echo "Pick a project: nix develop .#73mGe | .#78mBr | .#YAG-PreProcessing | .#YAP-PSD"
            '';
          };
        };
      }
    );
}
