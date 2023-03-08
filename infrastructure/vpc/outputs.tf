output "security_groups" {
  value = [aws_security_group.security-group.id, aws_security_group.internal-security-group.id]
}

output "subnet_id" {
  value = aws_subnet.public-subnet.id
}

output "key_name" {
  value = aws_key_pair.kp.key_name
}
