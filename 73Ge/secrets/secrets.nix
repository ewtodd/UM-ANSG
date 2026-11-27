let
  ethan-desktop = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIDlbs+h9OqZMIAC6b3i4tUcXC4PidfBFEQNdwrLS8g9G";
  ethan-laptop = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIF5aPPhXY+RssvL9znCFwHjkmUdi4KQkNSnAgd+AQqqx";
  allKeys = [
    ethan-desktop
    ethan-laptop
  ];
in
{
  "BEF.hpp.age".publicKeys = allKeys;
  "BEF.cpp.age".publicKeys = allKeys;
  "ConvertBEF.cpp.age".publicKeys = allKeys;
}
