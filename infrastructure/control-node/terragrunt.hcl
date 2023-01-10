terraform {
  source = "git::git@github.com:ntedeschiSG/ansible-modules.git//instance"
}

dependency "vpc" {
  config_path = "../vpc"
}

inputs = {
  instance_name   = "ntedeschi-postal-control"
  ami             = "ami-0885a26f897490207"
  instance_type   = "t3.medium"
  root_vol_size   = 30 
  key_name        = dependency.vpc.outputs.key_name

  security_groups = dependency.vpc.outputs.security_groups
  subnet_id       = dependency.vpc.outputs.subnet_id
  eip             = "eipalloc-0335c1dad94a79659"
  private_ip      = "10.0.1.0"
}
