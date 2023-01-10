terraform {
  source = "git::git@github.com:ntedeschiSG/terraform-modules.git//managed-node"
}

dependency "vpc" {
  config_path = "../../vpc"
}

inputs = {
  instance_name   = "ntedeschi-postal-1"
  #ami             = "ami-0c7c3448920e5e369" # NVIDIA GPU-Optimized AMI in region us-west-2
  ami             = "ami-003100cf7f2c7c335" # Rocky 9
  instance_type   = "m6i.2xlarge"
  root_vol_size   = 64 
  key_name        = dependency.vpc.outputs.key_name

  security_groups = dependency.vpc.outputs.security_groups
  subnet_id       = dependency.vpc.outputs.subnet_id
  private_ip      = "10.0.1.1"

  ebs_size        = 30 
  ebs_name        = "ntedeschi-postal-1"
  ebs_availability_zone = "us-west-2b"
}

