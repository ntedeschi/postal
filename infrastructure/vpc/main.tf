terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 4.34.0"
    }
  }

  required_version = ">= 1.3"
}

provider "aws" {
  region = var.region
}

resource "aws_key_pair" "kp" {
  key_name = var.key_name
  public_key = var.public_key
}

resource "aws_vpc" "vpc" {
  cidr_block = "10.0.0.0/16"
  enable_dns_hostnames = true
  enable_dns_support = true

  tags = {
    Name = var.vpc_name
  }
}

resource "aws_security_group" "security-group" {
  name = "postal-security-group"
  vpc_id = aws_vpc.vpc.id

  ingress {
    description = "ssh access"
    cidr_blocks =  var.cidr_blocks
    from_port = 22
    to_port = 22
    protocol = "tcp"
  }

  egress {
    from_port = 0
    to_port = 0
    protocol = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = {
    Name = var.vpc_name
  }
}

resource "aws_security_group" "internal-security-group" {
  name = "postal-internel-security-group"
  vpc_id = aws_vpc.vpc.id

  ingress {
    description = "ssh access"
    cidr_blocks = ["10.0.0.0/16"] 
    from_port = 22
    to_port = 22
    protocol = "tcp"
  }

  egress {
    from_port = 0
    to_port = 0
    protocol = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  tags = {
    Name = var.vpc_name
  }
}

resource "aws_internet_gateway" "ig" {
  vpc_id = aws_vpc.vpc.id
  tags = {
    Name = var.vpc_name
  }
}

resource "aws_route_table" "rt" {
  vpc_id = aws_vpc.vpc.id

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.ig.id
  }

  tags = {
    Name = var.vpc_name
  }
}

resource "aws_main_route_table_association" "rta" {
  vpc_id = aws_vpc.vpc.id
  route_table_id = aws_route_table.rt.id
}

resource "aws_subnet" "public-subnet" {
  vpc_id = aws_vpc.vpc.id
  cidr_block = "10.0.0.0/20"
  availability_zone = var.availability_zone
  map_public_ip_on_launch = true

  tags = {
    Name = var.vpc_name
  }
}
