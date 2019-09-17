add_executable(
	dipoleamplitude	
	main.cpp
	)


include_directories(../src)
target_link_libraries(
	dipoleamplitude
	PRIVATE
	mzipsat 
)


