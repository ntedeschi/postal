terraform {
  source = "git::git@github.com:ntedeschiSG/ansible-modules.git//instance"
}

dependency "vpc" {
  config_path = "../vpc"
}

inputs = {
  instance_name   = "postal-control"
  ami             = "ami-0672af4b5c29cece0" # Rocky 9  9.1.20230215
  instance_type   = "t3.medium"
  root_vol_size   = 30 
  key_name        = dependency.vpc.outputs.key_name

  security_groups = dependency.vpc.outputs.security_groups
  subnet_id       = dependency.vpc.outputs.subnet_id
  eip             = "eipalloc-062846e4c3093d37a"
  private_ip      = "10.0.1.0"
}
