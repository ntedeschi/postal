terraform {
  source = "git::git@github.com:ntedeschiSG/terraform-modules.git//managed-node"
}

dependency "vpc" {
  config_path = "../../vpc"
}

inputs = {
  instance_name   = "ntedeschi-postal-1"
  ami             = "ami-0672af4b5c29cece0" # Rocky 9  9.1.20230215
  instance_type   = "m6i.2xlarge"
  root_vol_size   = 64 
  key_name        = dependency.vpc.outputs.key_name

  security_groups = dependency.vpc.outputs.security_groups
  subnet_id       = dependency.vpc.outputs.subnet_id
  private_ip      = "10.0.1.1"

  ebs_size        = 128 
  ebs_name        = "postal-1"
  ebs_availability_zone = "us-west-2b"
}

