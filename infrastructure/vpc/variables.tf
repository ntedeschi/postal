variable "vpc_name" {
  default = "postal-vpc"
}

variable cidr_blocks {
  default = ["70.95.118.178/32", "91.219.212.195/32"]
}

variable "region" {
  description = "Oregon region"
  default     = "us-west-2"
}

variable "availability_zone" {
  default = "us-west-2b"  # 2b or 2c are required for p3 instances
}

variable "key_name" {
  default = "postal-key"
}

variable "public_key" {
  default = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAINBhsGV3MaUheDL/PzR1GY/QvWmk5gpqhNBAUdjq4sjq ntedeschi@protonmail.com"
}
