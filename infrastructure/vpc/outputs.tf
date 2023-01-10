output "security_groups" {
  value = [aws_security_group.sg-home.id, aws_security_group.sg-internal.id]
}

output "subnet_id" {
  value = aws_subnet.public-subnet.id
}

output "key_name" {
  value = aws_key_pair.kp.key_name
}
