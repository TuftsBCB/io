all:
	go install ./...

push:
	git push origin master
	git push github master
	git push tufts master

